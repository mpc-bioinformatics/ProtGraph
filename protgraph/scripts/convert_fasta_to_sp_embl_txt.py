import argparse
import os

import tqdm

from protgraph.cli import check_if_file_exists


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Small script to create a SP-EMBL-txt-Entries from FASTA-entries."
        "This can be very useful for FASTA-database and -entries, which are NOT in UniProt but should be utilized by ProtGraph."
        "Currently only a plain conversion from FASTA to SP-EMBL is provided. However this could be extended to include "
        "further feature information to peptides (like variants) using a csv entry. "
        "Note: The header should have 3 section '<pre>|<accession>|<description>', whre the accession is unique for the whole FASTA-file."
    )

    # Number of entries in fasta (for tqdm)
    parser.add_argument(
        "--num_entries", "-n", type=int, default=None,
        help="Number of entries in the fasta files. if provided, it can give an estimate an running time."
    )

    # Base Folder of generated Pickle Files
    parser.add_argument(
        "fasta_file", type=check_if_file_exists, nargs=1,
        help="Fasta file, where the header should be replaced"
    )

    # Feature-Table-Mapping
    # parser.add_argument(
    #     "--feature_tables", "-ft", type=str, default=None,
    #     help="Feature-Table-Mapping, to add additional information into the SP-EMBL, TODO"
    # )

    # Output SP-EMBL-file
    parser.add_argument(
        "--output", "-o", type=str, default="output.txt",
        help="Output fasta file. DEFAULT 'output.fasta' (NOTE: This file WILL be overwritten)"
    )

    return parser.parse_args()


def main():
    # Parse args
    args = parse_args()

    # Set in and output
    in_fasta = os.path.abspath(args.fasta_file[0])
    out_txt = os.path.abspath(args.output)
    # feature_mapping = os.path.abspath(args.feature_tables)

    # Parameters which can be set
    pass # PLACEHOLDER
