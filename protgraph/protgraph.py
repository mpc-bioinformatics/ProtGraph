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
    prot_graph(**ARGS)


def prot_graph(**kwargs):
    """ MAIN DESCRIPTION TODO """
    prot_graph_args = get_defaults_args()  # Get default values
    prot_graph_args.update(kwargs)

    # Check for required Arguments!
    if "files" not in kwargs or kwargs["files"] is None:
        raise TypeError("missing argument 'files'")

    # Set up queues
    entry_queue = Queue(10000)
    statistics_queue = Queue(10000)
    common_out_file_queue = Queue(10000)

    # Get the number of processes.
    number_of_procs = \
        cpu_count() - 1 if prot_graph_args["num_of_processes"] is None else prot_graph_args["num_of_processes"]
    prot_graph_args["num_of_processes"] = number_of_procs

    # Create processes
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

    # Start processes/threads in reverse!
    common_out_thread.start()
    main_write_thread.start()
    for p in graph_gen:
        p.start()
    entry_reader.start()

    # Check processes and Threads in Reverse
    graph_gen_stop_sent = False
    main_write_threads_stop_sent = False
    while True:
        time.sleep(1)

        # Are writing threads alive?
        if not main_write_thread.is_alive() and not common_out_thread.is_alive():
            # Then exit the program
            break

        # Are all consumers still alive?
        if all([not x.is_alive() for x in graph_gen]) and not main_write_threads_stop_sent:
            # Add None to the last queue to stop thread
            statistics_queue.put(None)
            common_out_file_queue.put(None)
            main_write_threads_stop_sent = True
            continue

        # Is producer still alive?
        if not entry_reader.is_alive() and not graph_gen_stop_sent:
            # Add None, to stop all processes
            for _ in range(number_of_procs):
                entry_queue.put(None)
            graph_gen_stop_sent = True
            continue


def format_help(parser, groups=None):
    """ This function only prints specific parts of the help message! """
    formatter = parser._get_formatter()

    formatter.add_usage(parser.usage, [y for x in groups for y in x._group_actions],
                        parser._mutually_exclusive_groups)

    # description
    formatter.add_text(parser.description)

    if groups is None:
        groups = parser._action_groups

    # positionals, optionals and user-defined groups
    for action_group in groups:
        formatter.start_section(action_group.title)
        formatter.add_text(action_group.description)
        formatter.add_arguments(action_group._group_actions)
        formatter.end_section()

    # epilog
    formatter.add_text(parser.epilog)

    # determine help from format above
    return formatter.format_help()


def create_parser():
    """ Creates the argument parser, using all methods from cli.py """
    # Set parser without help
    parser = argparse.ArgumentParser(
        description="ProtGraph: a graph generator for proteins and/or peptides and exporter to various formats",
        add_help=False
    )

    # Set custom HelpAction to execute (Callback)
    class HelpAction(argparse.Action):
        def __init__(self, option_strings, dest, nargs=0, **kwargs):
            super(HelpAction, self).__init__(option_strings, dest, nargs=0, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            print(format_help(parser, parser._action_groups[:3]))  # Here we shorten!
            parser.exit()

    # Add shortened help message
    parser.add_argument(
        "--help", "-h", action=HelpAction, default=argparse.SUPPRESS, help="Show the shortened help message"
    )
    cli.add_main_args(parser)

    # Get group names and the corresponding functions from cli.py
    cli_groups = [
        ("graph_generation", cli.add_graph_generation),
        ("statistics", cli.add_statistics),
        ("graph_exports", cli.add_graph_exports),
        ("redis_graph_export", cli.add_redis_graph_export),
        ("postgres_graph_export", cli.add_postgres_graph_export),
        ("mysql_graph_export", cli.add_mysql_graph_export),
        ("cassandra_graph_export", cli.add_cassandra_export),
        ("postgres_peptide_export", cli.add_postgres_peptide_export),
        ("mysql_peptide_export", cli.add_mysql_peptide_export),
        ("citus_peptide_export", cli.add_citus_peptide_export),
        ("fasta_peptide_export", cli.add_fasta_peptide_export),
        ("gremlin_graph_export", cli.add_gremlin_graph_export),
    ]

    # First add help Option for all parameters.
    helpgroup = parser.add_argument_group("Enter one of these flags for detailed information")
    helpgroup.add_argument(
        "--help_all", action='help', help="Show the complete help message for all possible arguments"
    )

    # Use another HelpAction to print the corresponding section (Callback)
    class DetailHelpAction(argparse.Action):
        def __init__(self, option_strings, dest, nargs=0, **kwargs):
            super(DetailHelpAction, self).__init__(option_strings, dest, nargs=0, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            # Get idx and print the first part of help WITH the corresponding help section
            idx = [x[0] for x in cli_groups].index(option_string[len("--help_"):])
            print(format_help(parser, [*parser._action_groups[:2], parser._action_groups[idx+3]]))
            parser.exit()

    # Set arguments in parser, including the HelpAction for each section
    for name, func in cli_groups:
        helpgroup.add_argument("--help_{}".format(name), action=DetailHelpAction, default=argparse.SUPPRESS)
        group = parser.add_argument_group(name)
        func(group)

    # Return the parser
    return parser


def parse_args(args=None):
    """ Get the parser and parse the args. Returning a dict of arguments """
    parser = create_parser()

    # Parse arguments
    args = parser.parse_args(args)

    # Retrieve all parameters in a dict
    parsed_args = dict()
    for action in parser._actions:
        if action.default is not argparse.SUPPRESS:
            parsed_args[action.dest] = getattr(args, action.dest)

    return parsed_args


def get_defaults_args():
    """ Get the parser and retrieve the defaults. Returning a dict of arguments """
    parser = create_parser()

    # Retrieve all parameters in a dict
    defaults = dict()
    for action in parser._actions:
        if action.default is not argparse.SUPPRESS:
            defaults[action.dest] = action.default

    return defaults


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
                "Number of Mutagens",
                "Number of Conflicts",
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
