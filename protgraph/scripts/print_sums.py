import argparse
import ast
import csv
import math
import os
import sys
from operator import add

import tqdm

from protgraph.graph_statistics import _add_lists

csv.field_size_limit(sys.maxsize)


def _check_if_file_exists(s: str):
    """ checks if a file exists. If not: raise Exception """
    # TODO copied from prot_graph.py
    if os.path.isfile(s):
        return s
    else:
        raise Exception("File '{}' does not exists".format(s))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Graph-Generator for Proteins/Peptides and Exporter to various formats"
    )

    # Statistics file
    parser.add_argument(
        "input_csv", type=_check_if_file_exists, nargs=1,
        help="File containing the statistics output from ProtGraph "
        "(E.G.: Generated via '-cnp', '-cnpm' or '-cnph', ...). "
        "All columns beginning with 'num' or 'list' can be used."
    )

    # Number of entries in csv
    parser.add_argument(
        "--num_entries", "-n", type=int, default=None,
        help="Number of entries in csv. If provided, an estimate of needed time can be reported. Defaults to none"
    )

    # Number of entries in csv
    parser.add_argument(
        "--column_index", "-cidx", type=int, default=2,
        help="The Index of the column of the graph-statistics file which should be summed up. "
        "Defaults to 2 (num_var_seq (isoforms)) (The column for counting can be different depending on the layout)"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Parameters
    statistics_file = args.input_csv[0]
    num_entries = args.num_entries

    # Static parameter, always 11th entry in statistics file
    column_index = args.column_index

    # Open al files and sort them accordingly
    with open(statistics_file, "r") as in_file:
        # Initialize CSV reader
        csv_in = csv.reader(in_file)

        # Skip Header
        headers = next(csv_in)

        # Set up summation method depending on type
        try:
            while True:
                first = next(csv_in)[column_index]
                if not first:
                    continue
                parsed_entry = ast.literal_eval(first)
                if type(parsed_entry) is int:
                    exe_func = add
                    summation = 0
                    break
                elif type(parsed_entry) is list:
                    exe_func = _add_lists
                    summation = []
                    break
                else:
                    raise Exception("Column '{}' (on index {}) cannot be summed, type is not list or int".format(
                        headers[column_index], column_index)
                    )
        except StopIteration:
            raise Exception("Column '{}' (on index {}) cannot be summed. All cells are empty!".format(
                    headers[column_index], column_index)
                )

        # Sum all entries depending on type
        try:
            for row in tqdm.tqdm(csv_in, unit="rows", total=num_entries):
                if row[column_index]:
                    summation = exe_func(
                        summation, ast.literal_eval(row[column_index])
                    )
        except Exception:
            raise Exception("Column '{}' (on index {}) contains different typed/corrupted entries".format(
                headers[column_index], column_index)
            )

        # Print results
        if type(summation) is int:
            summation = [summation]
        idx_len = str(int(math.log10(len(summation)))+1)
        entry_len = str(int(math.log10(max(summation)))+1)
        print("Results from column '{}':\n".format(headers[column_index]))
        print("Sum of each entry")
        for idx, entry in enumerate(summation):
            print(("{:>" + idx_len + "}:  {:>" + entry_len + "}").format(idx, entry))

        print("\n\n\nCummulative sum")
        k = 0
        for idx, entry in enumerate(summation):
            k += entry
            print(("{:>" + idx_len + "}:  {:>" + entry_len + "}").format(idx, k))
