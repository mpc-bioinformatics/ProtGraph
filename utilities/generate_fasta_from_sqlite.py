import os
import argparse
import zlib


import tqdm
import apsw

def check_if_file_exists(s: str):
    """ checks if a folder exists. If not: raise Exception """
    if os.path.isfile(s):
        return s
    else:
        raise Exception("Filed '{}' does not exists".format(s))


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="Tool to generate fasta files from a generated filesysteme trie by ProtGraph. "
        "Protgraph should be called prior with the flag '-epeptrie'"
    )

    # Base Folder of generated Pickle Files
    parser.add_argument(
        "sqlite_file", type=check_if_file_exists, nargs=1,
        help="Sqlite database file od exported peptides."
    )

    # Output fasta file
    parser.add_argument(
        "--output_file", "-o", type=str, default="sqlite_peptides.fasta",
        help="Output fasta file. DEFAULT 'sqlite_peptides.fasta' (NOTE: File WILL be overwritten)"
    )

    return parser.parse_args()

def export_with_text(conn, output, num_of_entries):
    """ TODO """
    cur = conn.cursor()
    cur.execute("SELECT peptide, meta from peptide_meta;")

    pep_id = 0
    for pep, meta in tqdm.tqdm(cur, total=num_of_entries, unit="entries"):

        # Write Output-Entry            
        output_fasta.write(
            ">pg|ID_" + str(pep_id) + "|" + meta + "\n" + '\n'.join(pep[i:i+60] for i in range(0, len(pep), 60)) + "\n"
        )
        # Increase id
        pep_id += 1

    cur.close()


def export_with_blob(conn, output, num_of_entries):
    """ TODO """

    cur = conn.cursor()
    cur.execute("SELECT peptide from peptide_meta;")

    pep_id = 0
    for rowid in tqdm.tqdm(range(1,num_of_entries+1), total=num_of_entries):
        
        
        bl = conn.blobopen("main", "peptide_meta", "meta", rowid, False)
        binar_meta = bl.read()
        bl.close()
        meta = ""
        stream = binar_meta
        while stream:
            dco = zlib.decompressobj()
            meta += "," + dco.decompress(stream).decode("ascii")
            stream = dco.unused_data
        meta = meta[1:]

        bl = conn.blobopen("main", "peptide_meta", "peptide", rowid, False)
        binar_pep = bl.read()
        bl.close()
        pep = zlib.decompress(binar_pep).decode("ascii")


        # Write Output-Entry            
        output_fasta.write(
            ">pg|ID_" + str(pep_id) + "|" + meta + "\n" + '\n'.join(pep[i:i+60] for i in range(0, len(pep), 60)) + "\n"
        )
        # Increase id
        pep_id += 1



    # Getting blob entries:
    # cur.execute("""
    # Select ROWID from peptides_meta where peptide = x'78da0b020000530053';
    # """
    # )
    # print("=========================================================")
    # a = cur.fetchone()

    # bl = self.conn.blobopen("main", "peptides_meta", "meta", a[0], 1)
    # binar = bl.read()
    # bl.close()

    # stream = binar
    # while stream:
    #     dco = zlib.decompressobj()
    #     print(dco.decompress(stream))
    #     stream = dco.unused_data


if __name__ == "__main__":
    # Parse args
    args = parse_args()

    abs_file = os.path.abspath(args.sqlite_file[0])

    conn = apsw.Connection(abs_file)
    cur = conn.cursor()

    cur.execute("PRAGMA table_info(peptide_meta);")
    results = cur.fetchall()

    if results[0][2] == "TEXT":
        export_method = export_with_text
    else:  # == BLOB
        export_method = export_with_blob

    cur.execute("SELECT MAX(rowid) from  peptide_meta;")
    num_of_entries = cur.fetchone()[0]
    cur.close()



    # Start Export for text
    with open(args.output_file, "w") as output_fasta:
        export_method(conn, output_fasta, num_of_entries)

    conn.close()



