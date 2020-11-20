
from tqdm import tqdm
import csv
from Bio import SwissProt

#TODO remove Exclude List!

def read_embl(path_to_empl: str, queue):
    exclude_prots = "/home/luxii/git/PROVAR/exclude.csv"
    with open(exclude_prots) as in_f:
        csv_reader = csv.reader(in_f)
        exclude_list = [x[0] for x in list(csv_reader)]

    s = SwissProt.parse(path_to_empl)
    for entry in tqdm(s, total=192814, mininterval=0.5, unit="proteins"): 
        if entry.accessions[0] in exclude_list: # TODO this removes all worst cases at least for VARIANT !!!
            continue

        queue.put(entry)
