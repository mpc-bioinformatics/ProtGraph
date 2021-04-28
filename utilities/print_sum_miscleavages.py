import argparse
import ast
import csv
import math
import os
import sys
from operator import add

import tqdm

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
        help="File containing the statistics output from ProtGraph (generated via '-cnpm')"
    )

    # Number of entries in csv
    parser.add_argument(
        "--num_entries", "-n", type=int, default=None,
        help="Number of entries in csv. If provided, an estimate of needed time can be reported. Default None"
    )

    return parser.parse_args()


def add_lists(list_a, list_b):
    """ Appends 0 if list A or B is too short. Does the operation '+' to two lists (vectors, etc.) """
    if len(list_a) > len(list_b):
        # Swap elements if the other is larger
        t = list_a
        list_a = list_b
        list_b = t
    return list(map(
        add,
        list_a + [0]*(len(list_b) - len(list_a)),
        list_b
    ))


if __name__ == "__main__":
    args = parse_args()

    # Parameters
    statistics_file = args.input_csv[0]
    num_entries = args.num_entries

    # Static parameter, always 11th entry in statistics file
    mis_entry = 10

    # Open al files and sort them accordingly
    with open(statistics_file, "r") as in_file:
        # Initialize CSV reader
        csv_in = csv.reader(in_file)

        # Skip Header
        k = next(csv_in)

        # Sum all entries
        summation = []
        for row in tqdm.tqdm(csv_in, unit="rows", total=num_entries):
            in_list = ast.literal_eval(row[mis_entry])
            summation = add_lists(summation, in_list)

        # print results
        idx_len = str(int(math.log10(len(summation)))+1)
        entry_len = str(int(math.log10(max(summation)))+1)
        print("Sum of each entry")
        for idx, entry in enumerate(summation):
            print(("{:>" + idx_len + "}:  {:>" + entry_len + "}").format(idx, entry))

        print("\n\n\nCummulative Sum")
        k = 0
        for idx, entry in enumerate(summation):
            k += entry
            print(("{:>" + idx_len + "}:  {:>" + entry_len + "}").format(idx, k))
