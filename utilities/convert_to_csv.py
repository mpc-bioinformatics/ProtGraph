import argparse
import csv
import os
import sys

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
        "input_stdout", type=_check_if_file_exists, nargs=1,
        help="File containing the statistics output from ProtGraph (generated via '-cnph')"
    )

    # Number of entries in csv
    parser.add_argument(
        "--output_csv", "-o", type=str, default="output.csv",
        help="Output csv file"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Parameters
    stdout = args.input_stdout[0]
    out_csv_file = args.output_csv

    # Open al files and sort them accordingly
    with open(out_csv_file, "w") as out_file, open(stdout, "r") as in_file:
        # Initialize CSV writer
        csv_out = csv.writer(out_file)

        # Write header
        csv_out.writerow([
            "query", "protein", "method", "max_vars", "results", "time",
        ])

        # Initial values
        query = None

        # Sum all entries
        for row in tqdm.tqdm(in_file, unit="lines"):
            if row.startswith("Query: "):
                query = row[len("Query: "):-1]
            elif row.startswith("Protein: "):
                entries = row.split(",")
                csv_out.writerow([
                    query,
                    entries[0][len("Protein: "):],
                    entries[3][len(" method: "):] if len(entries) == 5 else entries[3][len(" method: "):-1],
                    entries[4][len(" max_vars"):-1] if len(entries) == 5 else None,
                    entries[2][len(" Results: "):],
                    entries[1][len(" Time: "):],
                ])
