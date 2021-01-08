
from tqdm import tqdm
import csv
from Bio import SwissProt


def read_embl(path_to_embls: list, num_of_entries: int, exclude_csv: str, queue):
    ''' Reads entries from a list of existing embl files '''
    if exclude_csv is None:
        # If no exclude csv is provided, we execute the reading without an if checking! (performance)
        for input_f in path_to_embls:
            # For each entry: try to read it and
            # add it to the queue
            try: 
                entries = SwissProt.parse(input_f)
                for entry in entries:
                    queue.put(entry)
            except Exception as e: 
                print("File '{}' could not be parsed and was excluded. Reason: {}".format(input_f, e))

    else:         
        # If a exclude csv is provided, then a simple if check is added (reduced performance)
        with open(exclude_csv) as in_f:
            # Read the contents of the csv
            csv_reader = csv.reader(in_f)
            exclude_list = [x[0] for x in list(csv_reader)]

            for input_f in path_to_embls:
                # For each entry: try to read it and
                # add it to the queue
                try: 
                    entries = SwissProt.parse(input_f)
                    for entry in entries:
                        if entry.accessions[0] in exclude_list: 
                            continue # This effectively skips an entry at the cost to check whether to skip in EACH entry!
                        queue.put(entry)
                except Exception as e: 
                    print("File '{}' could not be parsed and was excluded. Reason: {}".format(input_f, e))

