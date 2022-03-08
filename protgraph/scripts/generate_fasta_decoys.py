import argparse
import os
import random

import tqdm

from protgraph.cli import check_if_file_exists


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Small script to generate decoys in FASTA files"
        "This may be usefull for some generated fastas by protgraph, since the original headers "
        "may be to large or foreign for other decoy generators (since some always assume >sp or >tr at the beginning.."
    )

    # Number of entries in fasta (for tqdm)
    parser.add_argument(
        "--num_entries", "-n", type=int, default=None,
        help="Number of entries in the fasta files. if provided, it can give an estimate an running time."
    )

    # Base Folder of generated Pickle Files
    parser.add_argument(
        "fasta_file", type=check_if_file_exists, nargs=1,
        help="Fasta file, where decoys should be generated"
    )

    # Decoy String
    parser.add_argument(
        "--decoy_string", "-ds", type=str, default="DECOY_",
        help="Choose what decoy string should be added "
    )

    # Decoy Method
    parser.add_argument(
        "--decoy_method", "-dm", type=str, default="reverse", choices=["reverse", "shuffle"],
        help="Select the decoy method. Defaults to reverse (reversing the complete sequence)"
    )

    # Header Template
    parser.add_argument(
        "--header_template", "-t", type=str, default="{orig_pre}|{decoy_string}{orig_id}|{orig_desc}",
        help="Set the header template specifically. You can use: choose from the following variables: "
        "{orig_pre}, {orig_id}, {orig_desc} (from source fasta, splitted by '|' and {index} from target "
        "headers, indicating the index of the written header. NOTE for non-decoys the {decoy_string} defaults to ''"
        "DEFAULT: {orig_pre}|{decoy_string}{orig_id}|{orig_desc}"
    )

    # Output fasta file
    parser.add_argument(
        "--output_fasta", "-of", type=str, default="decoy.fasta",
        help="Output fasta file. DEFAULT 'decoy.fasta' (NOTE: This file WILL be overwritten)"
    )

    return parser.parse_args()


def reverse(sequence):
    return sequence[::-1]


def shuffle(sequence):
    return ''.join(random.sample(sequence, len(sequence)))


def main():
    # Parse args
    args = parse_args()

    # Set in and output
    in_fasta = os.path.abspath(args.fasta_file[0])
    out_fasta = os.path.abspath(args.output_fasta)

    if args.decoy_method == "reverse":
        decoy_method = reverse
    elif args.decoy_method == "shuffle":
        decoy_method = shuffle

    # Parameters which can be set
    template = args.header_template + "\n"
    index = 0

    with open(in_fasta, "r") as f_in_fasta, open(out_fasta, "w") as f_out_fasta:
        # Do for each entry
        pbar = tqdm.tqdm(total=args.num_entries, unit="entries")

        # Iterate via idx-1 and idx simultaneously
        header_line = next(f_in_fasta)
        for next_line in f_in_fasta:

            # Distinguish between entry and nonsense lines
            if header_line[:1] == ">":
                pbar.update()
                # Write header Info
                next_line_n_nl = header_line[:-1]

                # Get header
                orig_list = [""]*3
                splits = next_line_n_nl.split("|", 2)
                for idx, v in enumerate(splits):
                    orig_list[idx] = v

                # Get Sequence
                sequence = next_line[:-1]
                for seq_line in f_in_fasta:
                    if not seq_line.startswith(">"):
                        sequence += seq_line[:-1]
                    else:
                        header_line = seq_line
                        break

                # Generate Decoy string
                decoy_sequence = decoy_method(sequence)

                # Write original entry
                f_out_fasta.write(template.format(
                    index=index,
                    orig_pre=orig_list[0],
                    orig_id=orig_list[1],
                    orig_desc=orig_list[2],
                    decoy_string=""
                ))
                f_out_fasta.write(
                    '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60)) + "\n"
                )

                # Write DECOY entry
                f_out_fasta.write(template.format(
                    index=index,
                    orig_pre=orig_list[0],
                    orig_id=orig_list[1],
                    orig_desc=orig_list[2],
                    decoy_string=args.decoy_string
                ))
                f_out_fasta.write(
                    '\n'.join(decoy_sequence[i:i+60] for i in range(0, len(decoy_sequence), 60)) + "\n"
                )

                # Increase index
                index += 1
            else:
                header_line = next_line

        pbar.close()
