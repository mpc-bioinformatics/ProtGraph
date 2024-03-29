import argparse
import multiprocessing
import os
from functools import lru_cache
from multiprocessing import Process, cpu_count

import igraph
import psycopg
import tqdm

from protgraph.export.peptides.pep_fasta import PepFasta

PF = PepFasta()


def check_if_folder_exists(s: str):
    """ checks if a folder exists. If not: raise Exception """
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
        "base_export_folder", type=check_if_folder_exists, nargs=1,
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
    # Flag to export the compact form or all peptides
    parser.add_argument(
        "--compact", "-c", default=False, action="store_true",
        help="Set this flag to generate a more compact fasta file, generating entries having multiple proteins "
        "per peptide. This can reduce the size of the fasta by half and also reduce the overall runtime."
    )
    # Number of peptides to be processed at once
    parser.add_argument(
        "--batch_size", "-b", type=int, default=256,
        help="Number of entries which should be processed at once by a process"
    )
    # Number of processes to be used
    parser.add_argument(
        "--number_procs", "-np", type=int, default=None,
        help="Number of processes to be used during generation of the FASTA file"
    )
    return parser.parse_args()


def get_graph_file(base_dir, accession):
    """ gets the folder/file to the graph file, from the accession. """
    return os.path.join(
        base_dir,
        *[x for x in accession[:-1]],
        accession[-1] + ".pickle"  # TODO only pickle?, do we want to change this?
    )


@lru_cache(maxsize=25000)
def get_graph(base_dir, accession):
    """ Usage of lru_cache to maximize the speed on frequently called graphs. """
    return igraph.read(get_graph_file(base_dir, accession))


def convert_entry_to_fasta(peptide, part_headers, pep_id):
    # Generate header: >pg|ENTRY_ID|UNIPROT_ISO(START:END, mssclvg:MISS; QUALIFIERS)
    content = ">pg|ID_" + str(pep_id) + "|" + ",".join(part_headers)
    content += "\n" + '\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n"
    return content


def execute(in_q, out_q, base_folder, id_size, proc_id, num_procs):
    id_generator = PF.unique_id_gen(proc_id=proc_id, num_of_processes=num_procs)

    while True:
        rows = in_q.get()
        if rows is None:
            break
        # Export information per entry
        strings = []
        for row in rows:
            peptide, part_header = get_pep_and_header_def(row[0], row[1], base_folder)
            strings.append(
                convert_entry_to_fasta(peptide, [part_header], next(id_generator))
            )

        out_q.put("".join(strings))


def get_pep_and_header_def(path, accession, _base_folder):
    graph = get_graph(_base_folder, accession)
    peptide = "".join(graph.vs[path[1:-1]]["aminoacid"])
    edges = graph.get_eids(path=path)
    if "cleaved" in graph.es[edges[0]].attributes():
        misses = sum(filter(None, graph.es[edges]["cleaved"]))
    else:
        misses = -1

    acc = PF._get_accession_or_isoform(graph.vs[path[1]])
    start_pos = PF._get_position_or_isoform_position(graph.vs[path[1]])
    end_pos = PF._get_position_or_isoform_position(graph.vs[path[-2]], end=True)
    l_str_qualifiers = PF._get_qualifiers(graph, edges)
    quali_entries = ",".join(l_str_qualifiers)

    part_header = "".join(
        [
            acc, "(", str(start_pos), ":", str(end_pos), ",",
            "mssclvg:", str(misses),
            quali_entries,  ")"
        ]
    )

    return peptide, part_header


def execute_compact(in_q, out_q, base_folder, id_size, proc_id, num_procs):
    id_generator = PF.unique_id_gen(proc_id=proc_id, num_of_processes=num_procs)
    while True:
        rows = in_q.get()
        if rows is None:
            break
        # Generate dict of peptides
        strings = []
        for row in rows:
            d = dict()
            for path, accession in zip(row[0], row[1]):
                peptide, part_header = get_pep_and_header_def(path, accession, base_folder)
                if peptide not in d:
                    d[peptide] = [[]]
                d[peptide][0].append(part_header)

            # write dict entries to fasta
            for key, val, in d.items():
                strings.append(
                    convert_entry_to_fasta(key, val[0], next(id_generator))
                )

        out_q.put("".join(strings))


def batch(iterable, n):
    """ Batch entries into equally sized lists """
    try:
        while True:
            b = []
            for _ in range(n):
                b.append(next(iterable))
            yield b
    except StopIteration:
        yield b


def write_thread(queue, output_file, total_entries, b_size):
    """ Write thread to write into the fasta file"""
    pbar = tqdm.tqdm(total=total_entries, mininterval=0.5, unit="peptides")
    with open(output_file, "w") as fasta_out:
        while True:
            lines = queue.get()
            if lines is None:
                break
            fasta_out.write(lines)
            pbar.update(min(b_size, total_entries - pbar.n))


def main():
    # Parse args
    args = parse_args()

    # Initialize connection
    with psycopg.connect(
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

        # Iterate over each result from queue in parallel
        print("Waiting for database results...")

        in_queue = multiprocessing.Queue(1000)
        out_queue = multiprocessing.Queue(1000)

        number_of_procs = \
            cpu_count() - 1 if args.number_procs is None else args.number_procs
        pep_gen = [
            Process(
                target=exe, args=(in_queue, out_queue, args.base_export_folder[0], 1000000, x, number_of_procs)
            )
            for x in range(number_of_procs)
        ]
        for p in pep_gen:
            p.start()

        # Write output via extra thread
        main_write_thread = Process(
            target=write_thread,
            args=(out_queue, args.output_file, args.num_entries, args.batch_size,)
        )
        main_write_thread.start()

        # be generator as main thread
        with conn.cursor(binary=True, name="server_side_cursor_protgraph") as cursor:
            cursor.execute(query)
            while True:
                result = cursor.fetchmany(args.batch_size)
                if len(result) == 0:
                    break
                in_queue.put(
                    result
                )

        # Inform threads/processes to stop
        for _ in range(number_of_procs):
            in_queue.put(None)
        for p in pep_gen:
            p.join()

        out_queue.put(None)
        main_write_thread.join()


if __name__ == "__main__":
    main()
