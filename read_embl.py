
from tqdm import tqdm
import csv
from Bio import SwissProt

#TODO remove Exclude List!

def read_embl(path_to_empls: list, totals: list, queue):
    for input_file, total in zip(path_to_empls, totals):
        exclude_prots = "exclude.csv" # TODO REMOVE this
        with open(exclude_prots) as in_f:
            csv_reader = csv.reader(in_f)
            exclude_list = [x[0] for x in list(csv_reader)]

        s = SwissProt.parse(input_file)
        for entry in tqdm(s, total=total, mininterval=0.5, unit="proteins"):   
            if entry.accessions[0] in exclude_list: # TODO this removes all worst cases at least for VARIANT !!!
                continue

            queue.put(entry)
