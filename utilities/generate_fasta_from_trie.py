import os
import argparse

def check_if_folder_exists(s: str):
    """ checks if a folder exists. If not: raise Exception """
    if os.path.isdir(s):
        return s
    else:
        raise Exception("Folder '{}' does not exists".format(s))


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Tool to generate fasta files from a generated filesysteme trie by ProtGraph. "
        "Protgraph should be called prior with the flag '-epeptrie'"
    )

    # Base Folder of generated Pickle Files
    parser.add_argument(
        "base_export_folder", type=check_if_folder_exists, nargs=1,
        help="Base folder of the exported graphs in pickle (exported by the flag '-epeptrie')"
    )

    # Output fasta file
    parser.add_argument(
        "--output_file", "-o", type=str, default="exported_trie_peptides.fasta",
        help="Output fasta file. DEFAULT 'exported_trie_peptides.fasta' (NOTE: File WILL be overwritten)"
    )

    return parser.parse_args()


if __name__ == "__main__":
    # Parse args
    args = parse_args()


    abs_path = os.path.abspath(args.base_export_folder[0])
    len_abs_path = len(abs_path)
    pep_id = 0

    with open(args.output_file, "w") as output_fasta:

        # traverse root directory, and list directories as dirs and files as files
        for root, dirs, files in os.walk(abs_path):
            if len(files) >= 1:
                with open(os.path.join(root, files[0]), "r") as h_file:
                    header = h_file.read()[1:]
                peptide = root[len_abs_path:].replace(os.sep, "")

                output_fasta.write(">pg|ID_" + str(pep_id) + "|" + header + "\n")
                pep_id += 1
                output_fasta.write('\n'.join(peptide[i:i+60] for i in range(0, len(peptide), 60)) + "\n")
