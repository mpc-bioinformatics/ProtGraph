import argparse
import csv
import os
import time
from multiprocessing import Process, Queue, cpu_count
from threading import Thread

from tqdm import tqdm

import protgraph.cli as cli
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
    prot_graph_args["num_of_processes"] = number_of_procs

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
            target=generate_graph_consumer, args=(entry_queue, statistics_queue, common_out_file_queue, i),
            kwargs=prot_graph_args
        )
        for i in range(number_of_procs)
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


def format_help(self, groups=None):
    # self == parser
    formatter = self._get_formatter()

    formatter.add_usage(self.usage, [y for x in  groups for y in x._group_actions],
                        self._mutually_exclusive_groups)

    # description
    formatter.add_text(self.description)

    if groups is None:
        groups = self._action_groups

    # positionals, optionals and user-defined groups
    for action_group in groups:
        formatter.start_section(action_group.title)
        formatter.add_text(action_group.description)
        formatter.add_arguments(action_group._group_actions)
        formatter.end_section()

    # epilog
    formatter.add_text(self.epilog)

    # determine help from format above
    return formatter.format_help()


def parse_args(args=None):
    # Arguments for Parser
    parser = argparse.ArgumentParser(
        description="ProtGraph: a graph generator for proteins and/or peptides and exporter to various formats",
        add_help=False
    )
    

    class HelpAction(argparse.Action):
        def __init__(self, option_strings, dest, nargs=0, **kwargs):
            super(HelpAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            print(format_help(parser, parser._action_groups[:3]))
            parser.exit()

    # Add basic cli options
    parser.add_argument("--help", "-h", action=HelpAction, help="Show the shortened help message")
    cli.add_main_args(parser)



    # Get group names and the corresponding function
    cli_groups = [
        ("graph_generation", cli.add_graph_generation),
        ("statistics", cli.add_statistics),
        ("graph_exports", cli.add_graph_exports),
        ("redis_graph_export", cli.add_redis_graph_export),
        ("postgres_graph_export", cli.add_postgres_graph_export),
        ("mysql_graph_export", cli.add_mysql_graph_export),
        ("postgres_peptide_export", cli.add_postgres_peptide_export),
        ("mysql_peptide_export", cli.add_mysql_peptide_export),
        ("fasta_peptide_export", cli.add_fasta_peptide_export),
        ("citus_peptide_export", cli.add_citus_peptide_export),
        ("gremlin_graph_export", cli.add_gremlin_graph_export),
    ]
    

    helpgroup = parser.add_argument_group("Enter one of these flags for detailed information")
    helpgroup.add_argument("--help_all", action='help', help="Show the complete help message for all possible arguments")


    class DetailHelpAction(argparse.Action):
        def __init__(self, option_strings, dest, nargs=0, **kwargs):
            super(DetailHelpAction, self).__init__(option_strings, dest, nargs=0, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            idx = [x[0] for x in cli_groups].index(option_string[len("--help_"):])
            print(format_help(parser, [*parser._action_groups[:2], parser._action_groups[idx+3]]))
            parser.exit()    

    check_helps = []
    for name, func in cli_groups:
        helpgroup.add_argument("--help_{}".format(name), action=DetailHelpAction)
        group = parser.add_argument_group(name)
        func(group)
        check_helps.append("help_{}".format(name))


    # Parse Arguments
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
        export_csv=args.export_csv,
        export_large_csv=args.export_large_csv,
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
    header_dict = dict()

    while True:
        # Wait and get next result entry
        entry = queue.get()

        # check if entry is None
        # If it is, we stop
        if entry is None:
            break

        if entry[2] and entry[0] in header_dict:
            continue

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

        if entry[2] and entry[0] not in header_dict:
            header_dict[entry[0]] = True

    # Close all opened files
    for _, val in out_dict.items():
        val.close()


if __name__ == "__main__":
    main()
