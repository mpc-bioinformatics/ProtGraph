import argparse
import csv
import os
import time
from multiprocessing import Process, Queue, cpu_count
from threading import Thread

from tqdm import tqdm

from graph_generator import generate_graph_consumer
from read_embl import read_embl


def main():
    """ MAIN DESCRIPTION TODO """

    graph_gen_args = parse_args()

    entry_queue = Queue(1000)  # TODO limit? or not limit
    statistics_queue = Queue()

    number_of_procs = \
        cpu_count() - 1 if graph_gen_args["num_of_processes"] is None else graph_gen_args["num_of_processes"]

    # Create Processes
    entry_reader = Process(
        target=read_embl,
        args=(
            graph_gen_args["files"], graph_gen_args["num_of_entries"],
            graph_gen_args["exclude_accessions"], entry_queue,
        )
    )
    graph_gen = [
        Process(target=generate_graph_consumer, args=(entry_queue, statistics_queue), kwargs=graph_gen_args)
        for _ in range(number_of_procs)
    ]
    main_write_thread = Thread(
        target=write_output_csv_thread,
        args=(statistics_queue, graph_gen_args["output_csv"], graph_gen_args["num_of_entries"],)
    )

    # Start Processes/threads in reverse!
    for p in graph_gen:
        p.start()
    entry_reader.start()
    main_write_thread.start()

    # Check Processes and Threads in Reverse
    graph_gen_stop_sent = False
    main_write_thread_stop_send = False
    while True:
        time.sleep(1)

        # Is the writing thread alive?
        if not main_write_thread.is_alive():
            # Then exit the program
            break

        # Are Consumers alive?
        if all([not x.is_alive() for x in graph_gen]) and not main_write_thread_stop_send:
            # Add None to the last queue to stop thread
            statistics_queue.put(None)
            main_write_thread_stop_send = True
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


def parse_args():
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
        choices=["trypsin", "skip"],
        help="Set the digestion method. Default is set to Trypsin."
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
        "--output_csv", "-o", default=os.path.join(os.path.dirname(__file__), "protein_graph_statistics.csv"),
        type=str,
        help="Set the output file, which will contain information about the ProteinGaph (in csv) NOTE: "
        "It will write to 'protein_graph_statistics.csv' and overwrite if such a file exists. Default is "
        "set to write it in this projects folder"
    )

    # Arguments for exporting
    parser.add_argument(
        "--export_output_folder", "-eo", default=os.path.join(os.path.dirname(__file__), "exported_graphs"), type=str,
        help="Set the output folder to specify the dirctory of exported graphs (dot, graphml, gml) NOTE: It will "
        "overwrite exisiting files. Default is set to export in this projects folder"
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
        "--export_redisgraph", "-eredisg", default=False, action="store_true",
        help="Set this flag to export to a redis-server having the module RedisGraph."
        "NOTE: Currently the YY_BUFFER_SIZE is set to 1MB max large queries, you "
        "may encounter errors while adding graphs. To resolve this, you may need to build RedisGraph yourself"
        # TODO check on https://github.com/RedisGraph/RedisGraph/issues/1486
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

    args = parser.parse_args()

    # Graph generation arguments in dict:
    graph_gen_args = dict(
        files=args.files,
        num_of_entries=args.num_of_entries,
        exclude_accessions=args.exclude_accessions,
        # num of parallel processes
        num_of_processes=args.num_of_processes,
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
        output_csv=args.output_csv,
        # Export files in folder
        export_output_folder=args.export_output_folder,
        export_dot=args.export_dot,
        export_graphml=args.export_graphml,
        export_gml=args.export_gml,
        # Export RedisGraph
        export_redisgraph=args.export_redisgraph,
        redisgraph_host=args.redisgraph_host,
        redisgraph_port=args.redisgraph_port,
        redisgraph_graph=args.redisgraph_graph,
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


if __name__ == "__main__":
    main()
    # TODO theses are the numbers of nodes
    # and edged for complete uniprot. Remove it from the code!
    # Nodes Count: 71031227
    # Edges Count: 77997502

    # Optimization is possible here!
    # We could concat Nodes together as long as there is only one  in and out degree
    # This optimization can happen before!!! doing the weighting (we can use the topological sort for this)
