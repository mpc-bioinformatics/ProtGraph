from multiprocessing import Queue, Process


from read_embl import read_embl
from graph_generator import generate_graph_consumer
from variation_generator import get_next_variant
# from database_writer import insert_to_database
import time


def main():
    # TODO set parameters for start


    
    INPUT_FILES = [] # List of UniProt EMBL Files (either .dat or .txt)
    INPUT_FILES_ENTRIES = [] # List of Entries in the INPUT_FILES



    QUEUE_SIZE = 10000 # How many elements inside a Queue


    # Graph Generation (what to Skip what not to Skip)



    # Output Writer(s?)



import argparse
def parse_args():
    parser = argparse.ArgumentParser(description="Graph-Generator for Proteins/Peptides and Exporter to various formats")

    ### Argument Parsing 


    ## Arguments for Graph Generation
    parser.add_argument(
        "--skip_isoforms", "-si", default=False, action="store_true",
        help="Set this flag to exclude isoforms 'VAR_SEQ' (and possible modification on them like variations, etc...) from the FeatureTable"
    )    
    parser.add_argument(
        "--skip_variants", "-sv", default=False, action="store_true",
        help="Set this flag to exclude 'VARIANT' from the FeatureTable"
    )
    parser.add_argument(
        "--skip_init_met", "-sm", default=False, action="store_true",
        help="Set this flag to exclude the skipping of the initiator methionine ('INIT_M' in FeatureTable) for proteins"
    )
    parser.add_argument(
        "--skip_signal", "-ss", default=False, action="store_true",
        help="Set this flag to exclude skipping the signal peptide ('SIGNAL' in FeatureTable) of specific proteins"
    )




    # Arguments for Graph Processing/Digestion
    parser.add_argument(
        "--digestion", "-d", type=str.lower, default="trypsin",
        choices=["trypsin", "skip"],
        help="Set the digestion method"
    )







    args = parser.parse_args()

    # Graph Generation arguments in dict:
    graph_gen_args = dict(
        skip_isoforms= args.skip_isoforms,
        skip_variants= args.skip_variants,
        skip_init_met= args.skip_init_met,
        skip_signal  = args.skip_signal

    )


    return None, graph_gen_args, None


if __name__ == "__main__":

    _, graph_gen_args, _ = parse_args()


    entry_queue = Queue(15000)
    graph_queue = Queue(15000)
    prot_variation_queue = Queue(15000)
    output_queue = Queue(15000)








    input_file = ["/home/luxdo/git/uniprot-filtered-organism Homo+sapiens+(Human)+[9606] .txt"] # TODO parameterize
    input_file_totals = [194237]

    input_file = ["/home/luxii/git/variants_generator/e_coli.dat"] # TODO parameterize
    input_file_totals = [9434]

    # input_file = ["/hdd1tb/uniprot-filtered-organism Homo+sapiens+(Human)+[9606] .txt"]
    # input_file_totals = [192814]

    # input_file = [
    #     "/hdd1tb/uniprot_sprot.dat",
    #     "/hdd1tb/uniprot_trembl.dat"
    # ]
    # input_file_totals = [
    #     563552,
    #     195104019
    # ]

    # Create Processes
    entry_reader = Process(target=read_embl, args=(input_file, input_file_totals, entry_queue,))
    graph_gen = [Process(target=generate_graph_consumer, args=(entry_queue, graph_queue), kwargs=graph_gen_args) for _ in range(2)]
    # prot_var_gen = [Process(target=get_next_variant, args=(graph_queue, prot_variation_queue,)) for _ in range(8)]
    # database_writer = [Process(target=insert_to_database, args=(prot_variation_queue, output_queue,)) for _ in range(16)]
    


    # Start Processes in reverse!
    # for p in database_writer:
    #     p.start()
    # for p in prot_var_gen:
    #     p.start()
    for p in graph_gen:
        p.start()
    entry_reader.start()






    t = []
    while True: 
        
        try:
            # time.sleep(3)
            # entry = prot_variation_queue.get(timeout=600)
            # x = prot_variation_queue.get()
            # t.append( entry)
            entry = graph_queue.get(timeout=3)
            # (prot, count) = prot_variation_queue.get()
            # t.append( entry)

        except Exception:

            pass




    # Nodes Count: 71031227
    # Edges Count: 77997502

    # Optimization is possible here!
    # We could concat Nodes together as long as there is only one  in and out degree
    # This optimization can happen before!!! doing the weighting (we can use the topological sort for this)


inputs = "/raid_storage/uniprot/uniprot_sprot.dat"
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: 
            break
        yield b





with open(inputs, "r",encoding="utf-8",errors='ignore') as f:
    print (sum(bl.count("//") for bl in blocks(f)))