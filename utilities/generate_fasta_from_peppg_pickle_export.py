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
    parser = argparse.ArgumentParser(
        description="Graph-Generator for Proteins/Peptides and Exporter to various formats"
    )

    # Connection to postgres
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

    # Statistics file
    parser.add_argument(
        "base_export_folder", type=_check_if_folder_exists, nargs=1,
        help="Base folder of the exported graphs in pickle (exported by the flags '-epickle', '-edirs')"
    )

    # Output folder
    parser.add_argument(
        "--output_file", "-o", type=str, default="exported_peptides.fasta",
        help="Output path for the lower and upper csv. (Default: exported.fasta"
    )

    # Number of entries in db, if not set we retrieve it via count()
    parser.add_argument(
        "--num_entries", "-n", type=int, default=None,
        help="Number of entries in the database. It will be retrieved via 'count(*)' if not provided"
    )

    return parser.parse_args()


def get_graph_file(base_dir, accession):
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
    return None


@lru_cache(maxsize=1000)
def get_graph(base_dir, accession):
    return igraph.read(get_graph_file(base_dir, accession))


if __name__ == "__main__":
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
            # Get number of entries
            if args.num_entries is None:
                print("Retrieving the number of rows...")
                cursor.execute("SELECT count(*) FROM peptides_meta;")
                n = cursor.fetchone()[0]
                args.num_entries = n

        print("Executing statement...")
        with conn.cursor(name="fetch_peptides_size") as cursor:
            cursor.itersize = 250000

            # query = """
            # SELECT path, accession, weight FROM peptides_meta
            # LEFT JOIN peptides ON peptides.id = peptides_meta.peptides_id
            # LEFT JOIN accessions ON accessions.id = peptides_meta.accession_id
            # ORDER BY accession  -- do we need them sorted? This may take alot of time! TODO DL
            # """

            # cursor.execute(query)

            # with open(args.output_file, "w") as fasta_out:
            #     print("Iterating over results...")
            #     # Do next iterations:
            #     for row in tqdm.tqdm(cursor, mininterval=0.5, total=args.num_entries, unit="peptides", initial=0):
            #         # Export information
            #         accession = row[1]
            #         graph = get_graph(args.base_export_folder[0], accession)
            #         peptide = "".join(graph.vs[row[0][1:-1]]["aminoacid"])
            #         qualifiers = get_qualifiers(graph, row[0])

            #         content = ">lcl|PEPTIDE_" + accession \
            #             + "|PATH=" + "->".join([str(i) for i in row[0]]) \
            #             + "|QUALIFIERS=" + ",".join(qualifiers)
            #         content += "\n" + '\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n"
            #         fasta_out.write(content)

            # TODO Allow the user to set a flag!
            query = """
            SELECT json_agg(peptides_meta.path), ARRAY_AGG(accessions.accession) FROM peptides
            JOIN peptides_meta ON peptides.id = peptides_meta.peptides_id
            JOIN accessions ON accessions.id = peptides_meta.accession_id
            GROUP BY peptides.id
            """

            cursor.execute(query)

            with open(args.output_file, "w") as fasta_out:
                print("Iterating over results...")
                # Do next iterations:
                for row in tqdm.tqdm(cursor, mininterval=0.5, total=args.num_entries, unit="peptides", initial=0):

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

                    # write to file
                    for key, val, in d.items():
                        content = ">lcl|ACCESSIONS=" + ",".join(val[0]) \
                            + "|QUALIFIERS=" + ",".join(";".join(x) if x is not None else "" for x in val[1])
                        content += "\n" + '\n'.join(key[i:i+60] for i in range(0, len(key), 60)) + "\n"
                        fasta_out.write(content)
