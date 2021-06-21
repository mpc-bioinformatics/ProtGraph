import argparse
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
        help="File containing the statistics output from ProtGraph (generated via '-cnph')"
    )

    # Number of entries in csv
    parser.add_argument(
        "--num_entries", "-n", type=int, default=None,
        help="Number of entries in csv. If provided, an estimate of needed time can be reported. Default None"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Parameters
    statistics_file = args.input_csv[0]
    num_entries = args.num_entries

    # Static parameter, always 11th entry in statistics file
    path_entry = 11

    # Open al files and sort them accordingly
    with open(statistics_file, "r") as in_file:
        # Initialize CSV reader
        csv_in = csv.reader(in_file)

        # Skip Header
        k = next(csv_in)

        # Sum all entries
        summation = 0
        for row in tqdm.tqdm(csv_in, unit="rows", total=num_entries):
            next_pos_paths = int(row[path_entry])
            summation += next_pos_paths

        # print results
        print("Sum of each entry:")
        print(summation)
