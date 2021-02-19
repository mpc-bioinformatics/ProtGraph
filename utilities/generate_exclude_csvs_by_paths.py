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
        "input_csv", type=_check_if_file_exists, nargs=1,
        help="File containing the statistics output from ProtGraph (generated via '-cnp')"
    )

    # Output folder
    parser.add_argument(
        "--output_folder", "-o", type=str, default=os.getcwd(),
        help="Output path for the lower and upper csv. (Default: current working directory)"
    )
    # threshould
    parser.add_argument(
        "--threhshold", "-t", type=int, default=10000000,
        help="Threshold to distinguish entries which are lower or higher. (Default 10 000 000)"
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
    base_path = args.output_folder
    statistics_file = args.input_csv[0]
    amount = args.threhshold
    num_entries = args.num_entries

    # Static parameter, always the 0th and 9th entry in statistics
    accession_entry = 0
    paths_entry = 9

    # Output files
    lower_file = os.path.join(base_path, "lower_by_paths.csv")
    upper_file = os.path.join(base_path, "upper_by_paths.csv")

    # Open al files and sort them accordingly
    with open(lower_file, "w") as l_in, open(upper_file, "w") as u_in, open(statistics_file, "r") as in_file:
        # Initialize CSV writer/reader
        l_csv_in, u_csv_in = csv.writer(l_in), csv.writer(u_in)
        csv_in = csv.reader(in_file)

        # Write header
        l_csv_in.writerow(["Accession"])
        u_csv_in.writerow(["Accession"])

        # Skip Header
        k = next(csv_in)

        paths_sum_lower = 0
        paths_sum_upper = 0
        prot_count_lower = 0
        prot_count_upper = 0
        for row in tqdm.tqdm(csv_in, unit="rows", total=num_entries):
            e = int(row[paths_entry])
            if e <= amount:
                l_csv_in.writerow([row[accession_entry]])
                paths_sum_lower += e
                prot_count_lower += 1
                continue

            u_csv_in.writerow([row[accession_entry]])
            paths_sum_upper += e
            prot_count_upper += 1

        print("Sum of paths which are lower  then {} (in total {} proteins):\n\t{}"
              .format(amount, prot_count_lower, paths_sum_lower))
        print("Sum of paths which are higher then {} (in total {} proteins):\n\t{}"
              .format(amount, prot_count_upper, paths_sum_upper))
