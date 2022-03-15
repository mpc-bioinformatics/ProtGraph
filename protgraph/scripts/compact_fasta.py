from collections import defaultdict
from tqdm import tqdm
import argparse
import os


from protgraph.cli import check_if_file_exists


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Small script to remove duplicated entries in a fasta and to summarize them into one"
        " This script works great in conjunction with the peptide fasta export from ProtGraph."
    )

    # Input Fasta
    parser.add_argument(
        "input_fasta", type=check_if_file_exists, nargs=1,
        help="Fasta file, the one which should be compacted"
    )

    # Custom Identifier
    parser.add_argument(
        "--identifier", "-id", type=str, default="ID_",
        help="Choose the identifier. Defaults to 'ID_'"
    )

    # Output fasta file
    parser.add_argument(
        "--output_fasta", "-o", type=str, default="compact.fasta",
        help="Output fasta file. DEFAULT 'compact.fasta' (NOTE: This file WILL be overwritten)"
    )

    return parser.parse_args()


def main():
    # Parse args
    args = parse_args()
    input_fasta = os.path.abspath(args.input_fasta[0])
    identifier = args.identifier
    output_fasta = os.path.abspath(args.output_fasta)

    # Read and Write files
    with open(input_fasta, "r") as inf, open(output_fasta, "w") as outf:
        header, sequence = None, ""
        seq_identifier = defaultdict(set)

        print("Generating dict of sets (to combine same sequences)...")
        for line in tqdm(inf, unit="lines"):

            if line.startswith(">"):
                if header is None:
                    # Edge Case very first entry 
                    header = line[:-1]
                    continue
                seq_identifier[sequence].add(header.split("|", 2)[2])
                sequence = ""
                header = line[:-1]    
            else:
                sequence += line[:-1]

        # Edge Case do it for the last entry
        seq_identifier[sequence].add(header.split("|", 2)[2])

        print("Writing new FASTA-file...")
        for index, (key, val) in tqdm(enumerate(seq_identifier.items()), unit="entries", total=len(seq_identifier)):
            outf.write(
                # Write header
                ">pg|" + identifier + str(index) + "|" + ",".join(val) + "\n" + 
                # Write sequence
                '\n'.join(key[i:i+60] for i in range(0, len(key), 60)) + "\n"
            )
