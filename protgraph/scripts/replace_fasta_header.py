import argparse
import os

import tqdm

from protgraph.cli import check_if_file_exists


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Small script to extract and replace fasta headers with a user-defined string. "
        "This may be usefull for some generated fastas by protgraph, since the original headers "
        "may be to large for e.g. search engines. Additionally the size of fasta could be reduced here!"
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

    # Output fasta file
    parser.add_argument(
        "--header_template", "-t", type=str, default="{orig_pre}|{index}|{orig_id}",
        help="Set the header template specifically. You can use: choose from the following variables: "
        "{orig_pre}, {orig_id}, {orig_desc} (from source fasta, splitted by '|' and {index} from target "
        "headers, indicating the index of the written header. "
        "DEFAULT: {orig_pre}|{index}|{orig_id}"
    )

    # Output fasta file
    parser.add_argument(
        "--output_fasta", "-of", type=str, default="output.fasta",
        help="Output fasta file. DEFAULT 'output.fasta' (NOTE: This file WILL be overwritten)"
    )

    parser.add_argument(
        "--output_header", "-oh", type=str, default="headers.txt",
        help="Output header file containing the headers from the original fasta. "
        "DEFAULT 'headers.txt' (NOTE: This file WILL be overwritten)"
    )

    return parser.parse_args()


def main():
    # Parse args
    args = parse_args()

    # Set in and output
    in_fasta = os.path.abspath(args.fasta_file[0])
    out_fasta = os.path.abspath(args.output_fasta)
    out_headers = os.path.abspath(args.output_header)

    # Parameters which can be set
    template = args.header_template + "\n"
    index = 0

    with open(in_fasta, "r") as f_in_fasta, \
         open(out_fasta, "w") as f_out_fasta, \
         open(out_headers, "w") as f_out_headers:
        # Do for each entry
        pbar = tqdm.tqdm(total=args.num_entries, unit="entries")
        for next_line in f_in_fasta:
            # Distinguish between header and entry
            if next_line[:1] == ">":
                pbar.update()
                # Write header Info
                next_line_n_nl = next_line[:-1]

                # Write new header
                orig_list = [""]*3
                splits = next_line_n_nl.split("|", 2)
                for idx, v in enumerate(splits):
                    orig_list[idx] = v
                f_out_fasta.write(template.format(
                    index=index,
                    orig_pre=orig_list[0],
                    orig_id=orig_list[1],
                    orig_desc=orig_list[2],
                ))

                # Extract and write original header
                f_out_headers.write(next_line)
                index += 1
            else:
                # Write entry info
                f_out_fasta.write(next_line)
        pbar.close()
