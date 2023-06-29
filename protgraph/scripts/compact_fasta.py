import argparse
import os
from collections import defaultdict

from tqdm import tqdm

from protgraph.cli import check_if_file_exists


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Small script to remove duplicated entries in fasta and to summarize them into one."
        " This script works great in conjunction with the peptide fasta export from ProtGraph and has been designed to"
        " work with huge FASTA-files (tested on ~200GB)"
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

    # Maximum number of entries in the dictionary
    parser.add_argument(
        "--entries_at_once", "-n", default=10_000_000, type=int,
        help="Number of entries to be checked at once. If set higher, more memory is consumed, "
        "but less time is needed. Default: 10 Million, which will need very roughly ~2-25GB memory "
        "(depending on the number of duplicates and the header sizes in the FASTA)."
    )

    return parser.parse_args()


def search_duplicates(
    dict_to_check: dict, in_file, cur_seek_position: int, initial_header: bytes, initial_header_pos: bytes
        ):
    """
    From the current position in file, look for duplicates for the remainder of this file.
    We update the dictionary where we check for duplicates. This method also returns the positions where
    a duplicate occured
    """
    # Initialize header and sequence
    header, sequence = initial_header, b""
    header_pos = initial_header_pos

    # Set for duplicated entries
    duplicated_entries = set()

    # Iterate over the remaining entries
    progress = tqdm(
        unit="B", total=in_file.seek(-1, os.SEEK_END) + 1 - cur_seek_position, position=1, leave=False, unit_scale=True
    )
    in_file.seek(cur_seek_position)
    line = in_file.readline()
    while line:
        progress.update(in_file.tell() - (progress.n + cur_seek_position))

        # On a new header, check if the previous entry is in the dict
        if line.startswith(b">"):
            if sequence in dict_to_check:
                # It is in the dict, include and update duplicated entries
                dict_to_check[sequence].add(header.split(b"|", 2)[2])
                duplicated_entries.add(header_pos)

            # Reset sequence and update header for next iteration
            sequence = b""
            header = line[:-1]
            header_pos = in_file.tell() - len(header)
        else:
            # Just append sequence
            sequence += line[:-1]

        # Next line for the next iteration
        line = in_file.readline()

    # For the last entry we check it seperately
    progress.update(in_file.tell() - (progress.n + cur_seek_position))
    if sequence in dict_to_check:
        # It is in the dict, include and update duplicated entries
        dict_to_check[sequence].add(header.split(b"|", 2)[2])
        duplicated_entries.add(header_pos)
    progress.close()

    # Return postions of duplicates
    return duplicated_entries


def write_dict(dict_to_write: dict, out_file, index_off=0, index_name="ID_"):
    """ Wrapper, to write  fasta entries of the entries in the dict """
    for index, (key, val) in tqdm(
        enumerate(dict_to_write.items()), unit="entries", total=len(dict_to_write), position=1, leave=False
    ):
        out_file.write(
            # Write header
            b">pg|" + index_name.encode() + str(index+index_off).encode() + b"|" + b",".join(val) + b"\n" +
            # Write sequence
            b'\n'.join(key[i:i+60] for i in range(0, len(key), 60)) + b"\n"
        )

    return index + index_off + 1  # plus one for the next iteration


def main():
    # Parse args
    args = parse_args()
    input_fasta = os.path.abspath(args.input_fasta[0])
    identifier = args.identifier
    output_fasta = os.path.abspath(args.output_fasta)
    entries_at_once = args.entries_at_once

    # Open Read and Write Files
    with open(input_fasta, "rb") as inf, open(output_fasta, "wb") as outf:
        # Set initial parameters
        header, sequence = None, b""
        seq_identifier = defaultdict(set)
        header_pos = None
        index_written = 0

        # Set set of entries of FASTA-entries to ignore
        ignore_duplicates_on = set()

        print("Searching for duplicates...")
        progress = tqdm(unit="B", total=inf.seek(-1, os.SEEK_END) + 1, position=0, unit_scale=True)

        # Iterate over each entry
        inf.seek(0)
        line = inf.readline()
        while line:
            progress.update(inf.tell() - progress.n)

            # If entry starts with header process the previous one
            if line.startswith(b">"):
                if header is None:
                    # Edge Case very first entry
                    header = line[:-1]
                    header_pos = inf.tell() - len(header)
                    line = inf.readline()
                    continue

                # Check if we ignore this entry
                if header_pos not in ignore_duplicates_on:
                    # If not, just add to the dictionary
                    seq_identifier[sequence].add(header.split(b"|", 2)[2])

                # Set sequence to none and set the new header (for the next iteration)
                sequence = b""
                header = line[:-1]
                header_pos = inf.tell() - len(header)

                # Check if the dict is already at its limit
                if len(seq_identifier) >= entries_at_once:
                    # If so, we go through the remainder of the file and write the dict

                    # Save position where we stopped to continue later
                    seek_to = inf.tell()  # Header is already saved

                    # Already delete duplicated "which we left behind" to save memory early on
                    ignore_duplicates_on = {x for x in ignore_duplicates_on if x > seek_to}

                    # Go through the remainder of the file to search for duplicates and save duplicated entries
                    ignore_duplicates_on.update(
                        search_duplicates(seq_identifier, inf, seek_to, header, header_pos)
                    )

                    # Go back to where we stopped
                    inf.seek(seek_to)

                    # Write entries in the dictionary
                    index_written = write_dict(seq_identifier, outf, index_written, index_name=identifier)

                    # Resest the dictionary for the next round of writing
                    seq_identifier.clear()
            else:
                # Just append the sequence
                sequence += line[:-1]

            # Read the next line for the next iteration
            line = inf.readline()

        # Edge Case: We simply include the last entry into the dict
        progress.update(inf.tell() - progress.n)
        progress.close()
        seq_identifier[sequence].add(header.split(b"|", 2)[2])

        print("Writing final FASTA-entries...")
        index_written = write_dict(seq_identifier, outf, index_written)
