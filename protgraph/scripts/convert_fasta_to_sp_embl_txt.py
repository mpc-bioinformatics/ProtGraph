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
        help="Number of entries in the fasta files. if provided, it can give an estimate of the running time."
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


def get_next_fasta_entry(fasta) -> tuple:
    """ Generator, returning parsed FASTA-entries """
    get_sequence = False # Flag to stop after the sequence was retrieved
    sequence = ""
    for line in fasta:
        # Iterate over each line of the FASTA-database
        if line.startswith(">"):
            # Case it is the header file
            if get_sequence:
                # We reached the next entry and can report the protein
                if sequence.isalpha() and sequence.isupper():
                    yield sequence, pre, accession, description
                    sequence = ""
                    get_sequence = False
                else:
                    print("WARNING: Entry {acc} has a malformed sequence".format(acc=accession))

            # Parse header information. Maybe we could extend this to regex?
            pre, accession, description = line[1:-1].split("|", 2)
            get_sequence = True

        else:
            # Simply append the sequences if we want to get it .
            if get_sequence:
                sequence += line[:-1]
    if sequence.isalpha() and sequence.isupper():
        yield sequence, pre, accession, description
    else:
        print("WARNING: Entry {acc} has a malformed sequence".format(acc=accession))


def generate_sp_embl_enty(man_sequence, man_accession, opt_feature):
    pass



def main():
    # Parse args
    args = parse_args()

    # Set in and output
    in_fasta = os.path.abspath(args.fasta_file[0])
    out_txt = os.path.abspath(args.output)
    # feature_mapping = os.path.abspath(args.feature_tables)


    with open(in_fasta, "r") as in_file:
        for sequence, pre, accession, description in get_next_fasta_entry(in_file):
            print(accession, sequence)
    

    # Parameters which can be set
    pass # PLACEHOLDER


if __name__ == "__main__":
    main()