import argparse
import csv
import os
import time
from multiprocessing import Process, Queue, cpu_count
from threading import Thread

from tqdm import tqdm

from protgraph.graph_generator import generate_graph_consumer
from protgraph.read_embl import read_embl


def main():
    ARGS = parse_args()
    prot_graph(ARGS)


def prot_graph(prot_graph_args):
    """ MAIN DESCRIPTION TODO """

    entry_queue = Queue(1000)  # TODO limit? or not limit
    statistics_queue = Queue()
    common_out_file_queue = Queue()

    number_of_procs = \
        cpu_count() - 1 if prot_graph_args["num_of_processes"] is None else prot_graph_args["num_of_processes"]

    # Create Processes
    entry_reader = Process(
        target=read_embl,
        args=(
            prot_graph_args["files"], prot_graph_args["num_of_entries"],
            prot_graph_args["exclude_accessions"], entry_queue
        )
    )
    graph_gen = [
        Process(
            target=generate_graph_consumer, args=(entry_queue, statistics_queue, common_out_file_queue),
            kwargs=prot_graph_args
        )
        for _ in range(number_of_procs)
    ]
    main_write_thread = Thread(
        target=write_output_csv_thread,
        args=(statistics_queue, prot_graph_args["output_csv"], prot_graph_args["num_of_entries"],)
    )
    common_out_thread = Thread(
        target=write_to_common_file,
        args=(common_out_file_queue,)
    )

    # Start Processes/threads in reverse!
    for p in graph_gen:
        p.start()
    entry_reader.start()
    main_write_thread.start()
    common_out_thread.start()

    # Check Processes and Threads in Reverse
    graph_gen_stop_sent = False
    main_write_threads_stop_send = False
    while True:
        time.sleep(1)

        # Is the writing thread alive?
        if not main_write_thread.is_alive() and not common_out_thread.is_alive():
            # Then exit the program
            break

        # Are Consumers alive?
        if all([not x.is_alive() for x in graph_gen]) and not main_write_threads_stop_send:
            # Add None to the last queue to stop thread
            statistics_queue.put(None)
            common_out_file_queue.put(None)
            main_write_threads_stop_send = True
            continue

        # Is Producer alive?
        if not entry_reader.is_alive() and not graph_gen_stop_sent:
            # Add None, to stop all processes
            for _ in range(number_of_procs):
                entry_queue.put(None)
            graph_gen_stop_sent = True
            continue


def check_if_file_exists(s: str):
    """ checks if a file exists. If not: raise Exception """
    if os.path.isfile(s):
        return s
    else:
        raise Exception("File '{}' does not exists".format(s))


def parse_args(args=None):
    # Arguments for Parser
    parser = argparse.ArgumentParser(
        description="Graph-Generator for Proteins/Peptides and Exporter to various formats"
    )

    # Needed Arguments for parsing (and additional information for it)
    parser.add_argument(
        "files", type=check_if_file_exists, nargs="+",
        help="Files containing the Swissprot/EMBL-Entries (either in .dat or .txt)"
    )
    parser.add_argument(
        "--num_of_entries", "-n", type=int, default=None,
        help="Number of entries across all files (summed). if given, it will an estimation of the runtime"
    )
    parser.add_argument(
        "--exclude_accessions", "-exclude", type=str, default=None,
        help="A csv file only containing accessions in the first row which should be excluded for processsing."
        " Setting this value may reduce the reading performance and therefore the throughput performance overall."
    )

    # Argument for number of Processes
    parser.add_argument(
        "--num_of_processes", "-np", type=int, default=None,
        help="The number of processes used to process entries. Each process can process an entry individually. "
        "The default value is 'cores - 1', since one process is reserved for reading entries"
    )

    # Flag to check if generated graphs are correctly generated
    parser.add_argument(
        "--verify_graph", "--verify", default=False, action="store_true",
        help="Set this flag to perform a check whether the graph was generated correctly. Here we explicitly check "
        "for parallel edges, for DAG and other properties."
    )

    # Arguments for graph generation
    parser.add_argument(
        "--skip_isoforms", "-si", default=False, action="store_true",
        help="Set this flag to exclude isoforms 'VAR_SEQ' (and possible modification on them like variations, etc...) "
        "from the FeatureTable"
    )
    parser.add_argument(
        "--skip_variants", "-sv", default=False, action="store_true",
        help="Set this flag to exclude 'VARIANT' from the FeatureTable"
    )
    parser.add_argument(
        "--skip_init_met", "-sm", default=False, action="store_true",
        help="Set this flag to exclude the skipping of the initiator methionine ('INIT_M' in "
        "FeatureTable) for proteins"
    )
    parser.add_argument(
        "--skip_signal", "-ss", default=False, action="store_true",
        help="Set this flag to exclude skipping the signal peptide ('SIGNAL' in FeatureTable) of specific proteins"
    )

    # Arguments for graph processing/digestion
    parser.add_argument(
        "--digestion", "-d", type=str.lower, default="trypsin",
        choices=["trypsin", "skip", "full"],
        help="Set the digestion method. The full digestion cleaves at every edge, which can be useful for retrieving "
        "all possible peptides with arbitrary cutting points. The digestion method skip skips the digestion "
        "completely. Default: Trypsin"
    )
    parser.add_argument(
        "--no_merge", "-nm", default=False, action="store_true",
        help="Set this flag to skip the merging process for chains of nodes and edges into a single node. Setting "
        "this option could drastically increase the size of the graph, especially its depth."
    )

    # Arguments for node and edge weights
    parser.add_argument(
        "--annotate_mono_weights", "-amw", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the monoisotopic weights. (Values are taken from "
        "the mass dictionary)"
    )
    parser.add_argument(
        "--annotate_avrg_weights", "-aaw", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the average weights. (Values are taken from "
        "the mass dictionary)"
    )
    parser.add_argument(
        "--annotate_mono_weight_to_end", "-amwe", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the monoisotopic end weights. This weight informs about "
        "how much weight is at least left to get to the end Node. NOTE: Applying this, also sets the monoisotopic "
        "weights"
    )
    parser.add_argument(
        "--annotate_avrg_weight_to_end", "-aawe", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the average end weights. This weight informs about "
        "how much weight is at least left to get to the end Node. NOTE: Applying this, also sets the average weights"
    )
    parser.add_argument(
        "--mass_dict_type", "-mdt",
        type=lambda s: int if s.lower() == "int" else (float if s.lower() == "float" else None), default="int",
        choices=[int, float], metavar="{int,float}",
        help="Set the type of the mass dictionary for amino acid. Default is set to int"
    )
    parser.add_argument(
        "--mass_dict_factor", "-mdf", type=float, default=1000000000,
        help="Set the factor for the masses inside the mass_dictionary. The default is set to 1 000 000 000, "
        "so that each mass can be converted into integers."
    )

    # Arguments for generation of graph statistics
    parser.add_argument(
        "--calc_num_possibilities", "-cnp", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This uses a dynamic programming approach to calculate this in an efficient manner."
    )
    parser.add_argument(
        "--calc_num_possibilities_miscleavages", "-cnpm", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This returns a list, sorted by the number of miscleavages (beginning at 0). "
        "Example: Returns: [1, 3, 5, 2] -> 1 path with 0 miscleavages, 3 paths with 1 miscleavage, 5 paths "
        "with 2 miscleavages, etc ... This uses a dynamic programming approach to calculate this in an efficient "
        "manner. NOTE: This may get memory heavy, depending on the proteins (especially on Titin)"
    )
    parser.add_argument(
        "--calc_num_possibilities_hops", "-cnph", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This returns a list, sorted by the number of hops (beginning at 0). "
        "Example: Returns: [0, 3, 5, 2] -> 0 paths with 0 hops, 3 paths with 1 hop, 5 paths "
        "with 2 hops, etc ... This uses a dynamic programming approach to calculate this in an efficient "
        "manner. NOTE: This mis even more memory heavy then binning on miscleavages. Of course it depends "
        "on the proteins (especially on Titin) NOTE: The dedicated start and end node is not counted here. "
        "If you traverse a graph, expect +2 more nodes in a path!"
    )

    parser.add_argument(
        "--output_csv", "-o", default=os.path.join(os.getcwd(), "protein_graph_statistics.csv"),
        type=str,
        help="Set the output file, which will contain information about the ProteinGaph (in csv) NOTE: "
        "It will write to 'protein_graph_statistics.csv' and overwrite if such a file exists. Default is "
        "set to the current working directory"
    )

    # Arguments for exporting
    parser.add_argument(
        "--export_output_folder", "-eo", default=os.path.join(os.getcwd(), "exported_graphs"), type=str,
        help="Set the output folder to specify the dirctory of exported graphs (dot, graphml, gml) NOTE: It will "
        "overwrite exisiting files. Default is set the current working directory"
    )
    parser.add_argument(
        "--export_in_directories", "-edirs", default=False, action="store_true",
        help="Set this flag to export files in directories (coded by accession) instead of directly by "
        "the accession name. This could be useful if millions of proteins are added into this tool, since "
        "a folder with millions of entries can be problematic in some cases."
    )
    parser.add_argument(
        "--export_dot", "-edot", default=False, action="store_true",
        help="Set this flag to export a dot file for each protein"
    )
    parser.add_argument(
        "--export_graphml", "-egraphml", default=False, action="store_true",
        help="Set this flag to export a GraphML file for each protein. This is the recommended export method."
    )
    parser.add_argument(
        "--export_gml", "-egml", default=False, action="store_true",
        help="Set this flag to export a GML file for each protein"
    )
    parser.add_argument(
        "--export_pickle", "-epickle", default=False, action="store_true",
        help="Set this flag to export a Pickle file for each protein"
    )
    parser.add_argument(
        "--export_redisgraph", "-eredisg", default=False, action="store_true",
        help="Set this flag to export to a redis-server having the module RedisGraph loaded."
    )
    parser.add_argument(
        "--redisgraph_host", type=str, default="localhost",
        help="Set the host name for the redis-server having the module RedisGraph. Default: localhost"
    )
    parser.add_argument(
        "--redisgraph_port", type=int, default=6379,
        help="Set the port for the redis-server having the module RedisGraph. Default: 6379"
    )
    parser.add_argument(
        "--redisgraph_graph", type=str, default="proteins",
        help="Set the graph name on the redis-server having the module RedisGraph. Default 'proteins'"
    )
    parser.add_argument(
        "--export_postgres", "-epg", default=False, action="store_true",
        help="Set this flag to export to a postgresql server."
        "NOTE: This will try to create the tables 'nodes' and 'edges' on a specified database."
        " Make sure the database in which the data should be saved also exists."
    )
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
        help="Set the database which will be used for the postgresql server. Default: proteins"
    )
    parser.add_argument(
        "--export_mysql", "-emysql", default=False, action="store_true",
        help="Set this flag to export to a MySQL server."
        "NOTE: This will try to create the tables 'nodes' and 'edges' on a specified database."
        " Make sure the database in which the data should be saved also exists."
    )
    parser.add_argument(
        "--mysql_host", type=str, default="127.0.0.1",
        help="Set the host name for the MySQL server. Default: 127.0.0.1"
    )
    parser.add_argument(
        "--mysql_port", type=int, default=3306,
        help="Set the port for the MySQL server. Default: 3306"
    )
    parser.add_argument(
        "--mysql_user", type=str, default="root",
        help="Set the username for the MySQL server. Default: root"
    )
    parser.add_argument(
        "--mysql_password", type=str, default="",
        help="Set the password for the MySQL server. Default: <empty>"
    )
    parser.add_argument(
        "--mysql_database", type=str, default="proteins",
        help="Set the database which will be used for the MySQL server. Default: proteins"
    )
    parser.add_argument(
        "--export_peptide_postgres", "-epeppg", default=False, action="store_true",
        help="Set this flag to export peptides (specifically paths) to a postgresql server."
        "NOTE: This will try to create the tables 'accessions' and 'peptides' on a specified database."
        " Make sure the database in which the data should be saved also exists. If problems occur, try "
        "to delete the generated tables and retry again."
    )
    parser.add_argument(
        "--pep_postgres_host", type=str, default="127.0.0.1",
        help="Set the host name for the postgresql server. Default: 127.0.0.1"
    )
    parser.add_argument(
        "--pep_postgres_port", type=int, default=5433,
        help="Set the port for the postgresql server. Default: 5433"
    )
    parser.add_argument(
        "--pep_postgres_user", type=str, default="postgres",
        help="Set the username for the postgresql server. Default: postgres"
    )
    parser.add_argument(
        "--pep_postgres_password", type=str, default="developer",
        help="Set the password for the postgresql server. Default: developer"
    )
    parser.add_argument(
        "--pep_postgres_database", type=str, default="proteins",
        help="Set the database which will be used for the postgresql server. Default: proteins"
    )
    parser.add_argument(
        "--pep_postgres_hops", type=int, default=None,
        help="Set the number of hops (max length of path) which should be used to get paths "
        "from a graph. NOTE: the larger the number the more memory may be needed. This depends on "
        "the protein which currently is processed. Default is set to 'None', so all lengths are considered."
    )
    parser.add_argument(
        "--pep_postgres_miscleavages", type=int, default=-1,
        help="Set this number to filter the generated paths by their miscleavages."
        "The protein graphs do contain infomration about 'infinite' miscleavages and therefor also return "
        "those paths/peptides. If setting (default) to '-1', all results are considered. However you can limit the "
        "number of miscleavages, if needed."
    )
    parser.add_argument(
        "--pep_postgres_skip_x",  default=False, action="store_true",
        help="Set this flag to skip to skip all entries, which contain an X"
    )
    parser.add_argument(
        "--pep_postgres_no_duplicates",  default=False, action="store_true",
        help="Set this flag to not insert duplicates into the database. "
        "NOTE: Setting this value decreases the performance drastically"
    )
    parser.add_argument(
        "--pep_postgres_use_igraph",  default=False, action="store_true",
        help="Set this flag to use igraph instead of netx. "
        "NOTE: If setting this flag, the peptide generation will be considerably faster "
        "but also consumes much more memory. Also, the igraph implementation DOES NOT go "
        "over each single edge, so some (repeating results) may never be discovered when using "
        "this flag."  # TODO see issue: https://github.com/igraph/python-igraph/issues/366
    )
    parser.add_argument(
        "--pep_postgres_min_pep_length", type=int, default=0,
        help="Set the minimum peptide length to filter out smaller existing path/peptides. "
        "Here, the actual number of aminoacid for a peptide is referenced. Default: 0"
    )
    parser.add_argument(
        "--pep_postgres_batch_size", type=int, default=25000,
        help="Set the batch size. This defines how many peptides are inserted at once. "
        "Default: 25000"
    )
    parser.add_argument(
        "--export_peptide_mysql", "-epepmysql", default=False, action="store_true",
        help="Set this flag to export peptides (specifically paths) to a MySQL server."
        "NOTE: This will try to create the tables 'accessions' and 'peptides' on a specified database."
        " Make sure the database in which the data should be saved also exists. If problems occur, try "
        "to delete the generated tables and retry again."
    )
    parser.add_argument(
        "--pep_mysql_host", type=str, default="127.0.0.1",
        help="Set the host name for the mysql server. Default: 127.0.0.1"
    )
    parser.add_argument(
        "--pep_mysql_port", type=int, default=3306,
        help="Set the port for the mysql server. Default: 3306"
    )
    parser.add_argument(
        "--pep_mysql_user", type=str, default="root",
        help="Set the username for the mysql server. Default: root"
    )
    parser.add_argument(
        "--pep_mysql_password", type=str, default="",
        help="Set the password for the mysql server. Default: ''"
    )
    parser.add_argument(
        "--pep_mysql_database", type=str, default="proteins",
        help="Set the database which will be used for the mysql server. Default: proteins"
    )
    parser.add_argument(
        "--pep_mysql_hops", type=int, default=None,
        help="Set the number of hops (max length of path) which should be used to get paths "
        "from a graph. NOTE: the larger the number the more memory may be needed. This depends on "
        "the protein which currently is processed. Default is set to 'None', so all lengths are considered."
    )
    parser.add_argument(
        "--pep_mysql_miscleavages", type=int, default=-1,
        help="Set this number to filter the generated paths by their miscleavages."
        "The protein graphs do contain infomration about 'infinite' miscleavages and therefor also return "
        "those paths/peptides. If setting (default) to '-1', all results are considered. However you can limit the "
        "number of miscleavages, if needed."
    )
    parser.add_argument(
        "--pep_mysql_skip_x",  default=False, action="store_true",
        help="Set this flag to skip to skip all entries, which contain an X"
    )
    parser.add_argument(
        "--pep_mysql_no_duplicates",  default=False, action="store_true",
        help="Set this flag to not insert duplicates into the database. "
        "NOTE: Setting this value decreases the performance drastically"
    )
    parser.add_argument(
        "--pep_mysql_use_igraph",  default=False, action="store_true",
        help="Set this flag to use igraph instead of netx. "
        "NOTE: If setting this flag, the peptide generation will be considerably faster "
        "but also consumes much more memory. Also, the igraph implementation DOES NOT go "
        "over each single edge, so some (repeating results) may never be discovered when using "
        "this flag."  # TODO see issue: https://github.com/igraph/python-igraph/issues/366
    )
    parser.add_argument(
        "--pep_mysql_min_pep_length", type=int, default=0,
        help="Set the minimum peptide length to filter out smaller existing path/peptides. "
        "Here, the actual number of aminoacid for a peptide is referenced. Default: 0"
    )
    parser.add_argument(
        "--pep_mysql_batch_size", type=int, default=25000,
        help="Set the batch size. This defines how many peptides are inserted at once. "
        "Default: 25000"
    )
    parser.add_argument(
        "--export_peptide_fasta", "-epepfasta", default=False, action="store_true",
        help="Set this flag to export peptides into a single fasta file."
    )
    parser.add_argument(
        "--pep_fasta_out", default=os.path.join(os.getcwd(), "peptides.fasta"),
        type=str,
        help="Set the output file for the peptide fasta export. Default: '${pwd}/peptides.fasta'. "
        "NOTE: This will overwrite existing files."
    )
    parser.add_argument(
        "--pep_fasta_hops", type=int, default=None,
        help="Set the number of hops (max length of path) which should be used to get paths "
        "from a graph. NOTE: the larger the number the more memory may be needed. This depends on "
        "the protein which currently is processed. Default is set to 'None', so all lengths are considered."
    )
    parser.add_argument(
        "--pep_fasta_miscleavages", type=int, default=-1,
        help="Set this number to filter the generated paths by their miscleavages."
        "The protein graphs do contain infomration about 'infinite' miscleavages and therefor also return "
        "those paths/peptides. If setting (default) to '-1', all results are considered. However you can limit the "
        "number of miscleavages, if needed."
    )
    parser.add_argument(
        "--pep_fasta_skip_x",  default=False, action="store_true",
        help="Set this flag to skip to skip all entries, which contain an X"
    )
    parser.add_argument(
        "--pep_fasta_use_igraph",  default=False, action="store_true",
        help="Set this flag to use igraph instead of netx. "
        "NOTE: If setting this flag, the peptide generation will be considerably faster "
        "but also consumes much more memory. Also, the igraph implementation DOES NOT go "
        "over each single edge, so some (repeating results) may never be discovered when using "
        "this flag."  # TODO see issue: https://github.com/igraph/python-igraph/issues/366
    )
    parser.add_argument(
        "--pep_fasta_min_pep_length", type=int, default=0,
        help="Set the minimum peptide length to filter out smaller existing path/peptides. "
        "Here, the actual number of aminoacid for a peptide is referenced. Default: 0"
    )
    parser.add_argument(
        "--pep_fasta_batch_size", type=int, default=25000,
        help="Set the batch size. This defines how many peptides are processed and written at once. "
        "Default: 25000"
    )
    parser.add_argument(
        "--export_peptide_citus", "-epepcit", default=False, action="store_true",
        help="Set this flag to export peptides (specifically paths) to a postgresql server."
        "NOTE: This will try to create the tables 'accessions' and 'peptides' on a specified database."
        " Make sure the database in which the data should be saved also exists. If problems occur, try "
        "to delete the generated tables and retry again."
    )
    parser.add_argument(
        "--pep_citus_host", type=str, default="127.0.0.1",
        help="Set the host name for the postgresql server with citus. Default: 127.0.0.1"
    )
    parser.add_argument(
        "--pep_citus_port", type=int, default=5433,
        help="Set the port for the postgresql server with citus. Default: 5433"
    )
    parser.add_argument(
        "--pep_citus_user", type=str, default="postgres",
        help="Set the username for the postgresql server with citus. Default: postgres"
    )
    parser.add_argument(
        "--pep_citus_password", type=str, default="developer",
        help="Set the password for the postgresql server with citus. Default: developer"
    )
    parser.add_argument(
        "--pep_citus_database", type=str, default="proteins",
        help="Set the database which will be used for the postgresql server with citus. Default: proteins"
    )
    parser.add_argument(
        "--pep_citus_hops", type=int, default=None,
        help="Set the number of hops (max length of path) which should be used to get paths "
        "from a graph. NOTE: the larger the number the more memory may be needed. This depends on "
        "the protein which currently is processed. Default is set to 'None', so all lengths are considered."
    )
    parser.add_argument(
        "--pep_citus_miscleavages", type=int, default=-1,
        help="Set this number to filter the generated paths by their miscleavages."
        "The protein graphs do contain infomration about 'infinite' miscleavages and therefor also return "
        "those paths/peptides. If setting (default) to '-1', all results are considered. However you can limit the "
        "number of miscleavages, if needed."
    )
    parser.add_argument(
        "--pep_citus_skip_x",  default=False, action="store_true",
        help="Set this flag to skip to skip all entries, which contain an X"
    )
    parser.add_argument(
        "--pep_citus_no_duplicates",  default=False, action="store_true",
        help="Set this flag to not insert duplicates into the database. "
        "NOTE: Setting this value decreases the performance drastically"
    )
    parser.add_argument(
        "--pep_citus_use_igraph",  default=False, action="store_true",
        help="Set this flag to use igraph instead of netx. "
        "NOTE: If setting this flag, the peptide generation will be considerably faster "
        "but also consumes much more memory. Also, the igraph implementation DOES NOT go "
        "over each single edge, so some (repeating results) may never be discovered when using "
        "this flag."  # TODO see issue: https://github.com/igraph/python-igraph/issues/366
    )
    parser.add_argument(
        "--pep_citus_min_pep_length", type=int, default=0,
        help="Set the minimum peptide length to filter out smaller existing path/peptides. "
        "Here, the actual number of aminoacid for a peptide is referenced. Default: 0"
    )
    parser.add_argument(
        "--pep_citus_batch_size", type=int, default=25000,
        help="Set the batch size. This defines how many peptides are inserted at once. "
        "Default: 25000"
    )
    parser.add_argument(
        "--export_gremlin", "-egremlin", default=False, action="store_true",
        help="Set this flag to export the graphs via gremlin to a gremlin server."
        "NOTE: The export is very slow, since it executes each node as a single query "
        "(tested on JanusGraph and Apache Gremlin Server). This exporter is not well implemented and may not work. "
        "This is due to difficulties implementing such an exporter in a global manner. "
        "To reduce the number of errors: Try to have a stable connection to the gremlin-server and also allocate "
        "enough resource for it, so that it can process the queries quick enough."
    )
    parser.add_argument(
        "--gremlin_url", type=str, default="ws://localhost:8182/gremlin",
        help="Set the url to the gremlin URL (no authentication). Default: 'ws://localhost:8182/gremlin'"
    )
    parser.add_argument(
        "--gremlin_traversal_source", type=str, default="g",
        help="Set the traversal source for remote. Default 'g'"
    )

    args = parser.parse_args(args)

    # Graph generation arguments in dict:
    graph_gen_args = dict(
        files=args.files,
        num_of_entries=args.num_of_entries,
        exclude_accessions=args.exclude_accessions,
        # num of parallel processes
        num_of_processes=args.num_of_processes,
        # verify graph flag
        verify_graph=args.verify_graph,
        # skip FTs?
        skip_isoforms=args.skip_isoforms,
        skip_variants=args.skip_variants,
        skip_init_met=args.skip_init_met,
        skip_signal=args.skip_signal,
        # Digestion and graph optiomization
        digestion=args.digestion,
        no_merge=args.no_merge,
        # Annotation arguments
        annotate_mono_weights=args.annotate_mono_weights,
        annotate_avrg_weights=args.annotate_avrg_weights,
        annotate_mono_weight_to_end=args.annotate_mono_weight_to_end,
        annotate_avrg_weight_to_end=args.annotate_avrg_weight_to_end,
        mass_dict_type=args.mass_dict_type,
        mass_dict_factor=args.mass_dict_factor,
        # Ouput CSV and num_of_paths arguments
        calc_num_possibilities=args.calc_num_possibilities,
        calc_num_possibilities_miscleavages=args.calc_num_possibilities_miscleavages,
        calc_num_possibilities_hops=args.calc_num_possibilities_hops,
        output_csv=args.output_csv,
        # Export files in folder
        export_output_folder=args.export_output_folder,
        export_in_directories=args.export_in_directories,
        export_dot=args.export_dot,
        export_graphml=args.export_graphml,
        export_gml=args.export_gml,
        export_pickle=args.export_pickle,
        # Export RedisGraph
        export_redisgraph=args.export_redisgraph,
        redisgraph_host=args.redisgraph_host,
        redisgraph_port=args.redisgraph_port,
        redisgraph_graph=args.redisgraph_graph,
        # Export postgresql
        export_postgres=args.export_postgres,
        postgres_host=args.postgres_host,
        postgres_port=args.postgres_port,
        postgres_user=args.postgres_user,
        postgres_password=args.postgres_password,
        postgres_database=args.postgres_database,
        # Export MySQL
        export_mysql=args.export_mysql,
        mysql_host=args.mysql_host,
        mysql_port=args.mysql_port,
        mysql_user=args.mysql_user,
        mysql_password=args.mysql_password,
        mysql_database=args.mysql_database,
        # Export peptides on postgresql
        export_peptide_postgres=args.export_peptide_postgres,
        pep_postgres_host=args.pep_postgres_host,
        pep_postgres_port=args.pep_postgres_port,
        pep_postgres_user=args.pep_postgres_user,
        pep_postgres_password=args.pep_postgres_password,
        pep_postgres_database=args.pep_postgres_database,
        pep_postgres_hops=args.pep_postgres_hops,
        pep_postgres_miscleavages=args.pep_postgres_miscleavages,
        pep_postgres_skip_x=args.pep_postgres_skip_x,
        pep_postgres_no_duplicates=args.pep_postgres_no_duplicates,
        pep_postgres_use_igraph=args.pep_postgres_use_igraph,
        pep_postgres_min_pep_length=args.pep_postgres_min_pep_length,
        pep_postgres_batch_size=args.pep_postgres_batch_size,
        # Export peptides on mysql
        export_peptide_mysql=args.export_peptide_mysql,
        pep_mysql_host=args.pep_mysql_host,
        pep_mysql_port=args.pep_mysql_port,
        pep_mysql_user=args.pep_mysql_user,
        pep_mysql_password=args.pep_mysql_password,
        pep_mysql_database=args.pep_mysql_database,
        pep_mysql_hops=args.pep_mysql_hops,
        pep_mysql_miscleavages=args.pep_mysql_miscleavages,
        pep_mysql_skip_x=args.pep_mysql_skip_x,
        pep_mysql_no_duplicates=args.pep_mysql_no_duplicates,
        pep_mysql_use_igraph=args.pep_mysql_use_igraph,
        pep_mysql_min_pep_length=args.pep_mysql_min_pep_length,
        pep_mysql_batch_size=args.pep_mysql_batch_size,
        # Export peptides to fasta file
        export_peptide_fasta=args.export_peptide_fasta,
        pep_fasta_out=args.pep_fasta_out,
        pep_fasta_hops=args.pep_fasta_hops,
        pep_fasta_miscleavages=args.pep_fasta_miscleavages,
        pep_fasta_skip_x=args.pep_fasta_skip_x,
        pep_fasta_use_igraph=args.pep_fasta_use_igraph,
        pep_fasta_min_pep_length=args.pep_fasta_min_pep_length,
        pep_fasta_batch_size=args.pep_fasta_batch_size,
        # Export peptides on postgresql
        export_peptide_citus=args.export_peptide_citus,
        pep_citus_host=args.pep_citus_host,
        pep_citus_port=args.pep_citus_port,
        pep_citus_user=args.pep_citus_user,
        pep_citus_password=args.pep_citus_password,
        pep_citus_database=args.pep_citus_database,
        pep_citus_hops=args.pep_citus_hops,
        pep_citus_miscleavages=args.pep_citus_miscleavages,
        pep_citus_skip_x=args.pep_citus_skip_x,
        pep_citus_no_duplicates=args.pep_citus_no_duplicates,
        pep_citus_use_igraph=args.pep_citus_use_igraph,
        pep_citus_min_pep_length=args.pep_citus_min_pep_length,
        pep_citus_batch_size=args.pep_citus_batch_size,
        # Export Gremlin (partially finished)
        export_gremlin=args.export_gremlin,
        gremlin_url=args.gremlin_url,
        gremlin_traversal_source=args.gremlin_traversal_source,
    )

    return graph_gen_args


def write_output_csv_thread(queue, out_file, total_num_entries):
    """
        The statistics writing thread, which writes to 'out_file', overwriting its
        contents if existing.
    """
    # Generate Progrssbar
    pbar = tqdm(total=total_num_entries, mininterval=0.5, unit="proteins")

    # (Over-)Write to out_file
    with open(out_file, "w") as out_f:
        csv_writer = csv.writer(out_f)

        # Write Header Row
        csv_writer.writerow(
            [
                "Accession",
                "Entry ID",
                "Number of isoforms",
                "Has INIT_MET",
                "Has SIGNAL",
                "Number of variants",
                "Number of cleaved edges",
                "Number of nodes",
                "Number of edges",
                "Num of possible paths",
                "Num of possible paths (by miscleavages 0, 1, ...)",
                "Num of possible paths (by hops 0, 1, ...)",
                "Protein description"
            ]
        )

        while True:
            # Wait and get next result entry
            next_entry = queue.get()

            # check if next_entry is None
            # If it is, we stop
            if next_entry is None:
                break

            # Write Data Row and update progress
            csv_writer.writerow(next_entry)
            pbar.update()

    # Close pbar, since we finished
    pbar.close()


def write_to_common_file(queue):
    """
    This method accepts a tuple from a queue where:
    the first element in tuple contains the destination of the file
    and the second the content which should be writton on
    """
    out_dict = dict()

    while True:
        # Wait and get next result entry
        entry = queue.get()

        # check if entry is None
        # If it is, we stop
        if entry is None:
            break

        # Write Data Row and update progress
        try:
            out_dict[entry[0]].write(entry[1])
        except Exception:
            # Create folder if needed
            if os.path.dirname(entry[0]) != "":
                os.makedirs(
                    os.path.dirname(entry[0]),
                    exist_ok=True
                )
            # Set entry
            out_dict[entry[0]] = open(entry[0], "w")
            # Rewrite first line!
            out_dict[entry[0]].write(entry[1])

    # Close all opened files
    for _, val in out_dict.items():
        val.close()


if __name__ == "__main__":
    main()
