import argparse
import os
from functools import lru_cache

import igraph
import psycopg2
import tqdm


def _check_if_folder_exists(s: str):
    """ checks if a file exists. If not: raise Exception """
    # TODO copied from prot_graph.py
    if os.path.isdir(s):
        return s
    else:
        raise Exception("Folder '{}' does not exists".format(s))


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Tool to generate fasta files from generated pickle graphs and peptide exports (Postgres). "
        "Protgraph should be called prior with the flags '-edirs', '-epickle' and'-epeppg'. "
    )

    # Parameters for the connection to postgres
    parser.add_argument(
        "--postgres_host", type=str, default="127.0.0.1",
        help="Set the host name for the postgresql server. Default: 127.0.0.1"
    )
    parser.add_argument(
        "--postgres_port", type=int, default=5433,
        help="Set the port for the postgresql server. Default: 5433"
    )
    parser.add_argument(
        "--postgres_user", type=str, default="postgres",
        help="Set the username for the postgresql server. Default: postgres"
    )
    parser.add_argument(
        "--postgres_password", type=str, default="developer",
        help="Set the password for the postgresql server. Default: developer"
    )
    parser.add_argument(
        "--postgres_database", type=str, default="proteins",
        help="Set the database which will be used on the postgresql server. Default: proteins"
    )

    # Base Folder of generated Pickle Files
    parser.add_argument(
        "base_export_folder", type=_check_if_folder_exists, nargs=1,
        help="Base folder of the exported graphs in pickle (exported by the flags '-epickle', '-edirs')"
    )

    # Output fasta file
    parser.add_argument(
        "--output_file", "-o", type=str, default="exported_peptides.fasta",
        help="Output fasta file. DEFAULT 'exported_peptides.fasta' (NOTE: File WILL be overwritten)"
    )

    # Number of entries in db, if not set we retrieve it via count()
    parser.add_argument(
        "--num_entries", "-n", type=int, default=None,
        help="Number of entries in the database. It will be retrieved via 'count(*)' if not provided"
    )
    # Flag to export the compact form or all peptide
    parser.add_argument(
        "--compact", "-c", default=False, action="store_true",
        help="Set this flag to generate a more compact fasta file, generating entries having multiple proteins "
        "per peptide. This can reduce the size of the fasta by half and also reduce the overall runtime."
    )
    return parser.parse_args()


def get_graph_file(base_dir, accession):
    """ gets the folder/file to the graph file, from the accession. """
    return os.path.join(
        base_dir,
        *[x for x in accession[:-1]],
        accession[-1] + ".pickle"  # TODO only pickle?, do we want to change this?
    )


def get_qualifiers(graph, path: list):
    """ Returns the qualifiers attributes (edges; in order, if existent) """
    # Get qualifier
    qualifiers = []
    for x, y in zip(path, path[1:]):
        edge = graph.es.find(_between=((x,), (y,)))
        # if it exits
        if "qualifiers" in edge.attributes():
            if edge["qualifiers"] is not None:
                for q in edge["qualifiers"]:
                    qualifiers.append(q.type)

    # return the aminoacids of the path
    if len(qualifiers) != 0:
        return qualifiers
    return None  # We return None if not present


@lru_cache(maxsize=25000)
def get_graph(base_dir, accession):
    """ Usage of lru_cache to maximize the speed on frequently called graphs. """
    return igraph.read(get_graph_file(base_dir, accession))


def write_entry_to_fasta(writer, peptide, accession_list, qualifier_list):
    content = ">lcl|ACCESSIONS=" + ",".join(accession_list) \
        + "|QUALIFIERS=" + ",".join(";".join(x) if x is not None else "" for x in qualifier_list)
    content += "\n" + '\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n"
    writer.write(content)


def execute(row, fasta):
    # Export information per entry
    accession = row[1]
    graph = get_graph(args.base_export_folder[0], accession)
    peptide = "".join(graph.vs[row[0][1:-1]]["aminoacid"])
    qualifiers = get_qualifiers(graph, row[0])
    write_entry_to_fasta(fasta, peptide, [accession], [qualifiers])


def execute_compact(row, fasta):
    # Generate dict of peptides
    d = dict()
    for path, accession in zip(row[0], row[1]):
        graph = get_graph(args.base_export_folder[0], accession)
        peptide = "".join(graph.vs[path[1:-1]]["aminoacid"])
        qualifiers = get_qualifiers(graph, path)
        if peptide not in d:
            d[peptide] = [[], []]
        d[peptide][0].append(accession)
        d[peptide][1].append(qualifiers)

    # write dict entries to fasta
    for key, val, in d.items():
        write_entry_to_fasta(fasta, key, val[0], val[1])


if __name__ == "__main__":
    # Parse args
    args = parse_args()

    # Initialize connection
    with psycopg2.connect(
            host=args.postgres_host,
            port=args.postgres_port,
            user=args.postgres_user,
            password=args.postgres_password,
            dbname=args.postgres_database
    ) as conn:
        with conn.cursor(name="fetch_peptides_size") as cursor:
            # Get number of entries if needed
            if args.num_entries is None:
                print("Retrieving the number of rows...")
                if args.compact:
                    cursor.execute("SELECT count(*) FROM peptides;")
                else:
                    cursor.execute("SELECT count(*) FROM peptides_meta;")
                n = cursor.fetchone()[0]
                args.num_entries = n

        print("Executing statement...")
        with conn.cursor(name="fetch_peptides_size") as cursor:
            # Set cursor iteration size (fetch x many entries at once)
            cursor.itersize = 250000

            # Get query function, depending on compact or not
            if args.compact:
                query = """
                SELECT json_agg(peptides_meta.path), ARRAY_AGG(accessions.accession) FROM peptides
                JOIN peptides_meta ON peptides.id = peptides_meta.peptides_id
                JOIN accessions ON accessions.id = peptides_meta.accession_id
                GROUP BY peptides.id
                """
                exe = execute_compact
            else:
                query = """
                SELECT path, accession FROM peptides_meta
                LEFT JOIN peptides ON peptides.id = peptides_meta.peptides_id
                LEFT JOIN accessions ON accessions.id = peptides_meta.accession_id
                """
                exe = execute

            # Execute query:
            cursor.execute(query)

            with open(args.output_file, "w") as fasta_out:
                print("Iterating over results...")
                # Execute for each row:
                for row in tqdm.tqdm(cursor, mininterval=0.5, total=args.num_entries, unit="peptides", initial=0):
                    exe(row, fasta_out)
