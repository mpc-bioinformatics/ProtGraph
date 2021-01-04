
from tqdm import tqdm
import csv
from Bio import SwissProt


def read_embl(path_to_embls: list, num_of_entries: int, queue):
    ''' Reads entries from a list of existing embl files '''
    for input_f in path_to_embls:
        # For each entry: try to read it and
        # add it to the queue
        try: 
            entries = SwissProt.parse(input_f)
            for entry in entries:
                queue.put(entry)
        except Exception as e: 
            print("File '{}' could not be parsed and was excluded. Reason: {}".format(input_f, e))

