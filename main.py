from multiprocessing import Queue, Process


from read_embl import read_embl
from graph_generator import generate_graph_consumer
from variation_generator import get_next_variant
# from database_writer import insert_to_database
import time

from tqdm import tqdm

def main():
    # TODO set parameters for start


    
    INPUT_FILES = [] # List of UniProt EMBL Files (either .dat or .txt)
    INPUT_FILES_ENTRIES = [] # List of Entries in the INPUT_FILES



    QUEUE_SIZE = 10000 # How many elements inside a Queue


    # Graph Generation (what to Skip what not to Skip)



    # Output Writer(s?)

def

import argparse
def parse_args():
    parser = argparse.ArgumentParser(description="Graph-Generator for Proteins/Peptides and Exporter to various formats")

    ### Argument parsing 
    parser.add_argument(
        "--output_csv", "-o", default="protein_graph_statistics.csv", type=argparse.types.Path(file=True, exists=False),
        help="Set the output file, which will contain information about the ProteinGaph (in csv)"
    )    

    ## Arguments for graph generation
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




    ## Arguments for graph processing/digestion
    parser.add_argument(
        "--digestion", "-d", type=str.lower, default="trypsin",
        choices=["trypsin", "skip"],
        help="Set the digestion method. Default is set to Trypsin."
    )
    parser.add_argument(
        "--no_merge", "-nm", default=False, action="store_true",
        help="Set this flag to skip the merging process for chains of nodes and edges "
        "into a single node. Setting this option could drastically increase the size of the graph, especially its depth."
    )


    ## Arguments for node and edge weights
    parser.add_argument(
        "--annotate_mono_weights", "-amw", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the monoisotopic weights. (Values are taken from the mass dictionary)"
    )
    parser.add_argument(
        "--annotate_avrg_weights", "-aaw", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the average weights. (Values are taken from the mass dictionary)"
    )
    parser.add_argument(
        "--annotate_mono_end_weights", "-amew", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the monoisotopic end weights. "
        "This weight informs about how much weight is at least left to get to the end Node. NOTE: Applying this, also sets the monoisotopic weights"
    )
    parser.add_argument(
        "--annotate_avrg_end_weights", "-aaew", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the average end weights. "
        "This weight informs about how much weight is at least left to get to the end Node. NOTE: Applying this, also sets the average weights"
    )
    parser.add_argument(
        "--mass_dict_type", "-mdt", type=lambda s: int if s.lower() == "int" else (float if s.lower() == "float" else None), default="int",
        choices=[int, float], metavar = "{int,float}",
        help="Set the type of the mass dictionary for amino acid. Default is set to int"
    )
    parser.add_argument(
        "--mass_dict_factor", "-mdf", type=float, default=1000000000,
        help="Set the factor for the masses inside the mass_dictionary. The default is set to 1 000 000 000, so that each mass can be converted into integers."
    )

    ## Arguments for generation of graph statistics
    parser.add_argument(
        "--calc_num_possibilities", "-cnp", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will be calculated. "
        "This uses a dynamic programming approach to calculate this in an efficient manner."
    )



    args = parser.parse_args()




    # Graph generation arguments in dict:
    graph_gen_args = dict(
        skip_isoforms= args.skip_isoforms,
        skip_variants= args.skip_variants,
        skip_init_met= args.skip_init_met,
        skip_signal  = args.skip_signal,

        digestion    = args.digestion,
        no_merge     = args.no_merge,

        annotate_mono_weights = args.annotate_mono_weights,
        annotate_avrg_weights = args.annotate_avrg_weights,
        annotate_mono_end_weights = args.annotate_mono_end_weights,
        annotate_avrg_end_weights = args.annotate_avrg_end_weights,
        mass_dict_type = args.mass_dict_type,
        mass_dict_factor = args.mass_dict_factor,

        calc_num_possibilities = args.calc_num_possibilities
    )


    return None, graph_gen_args, None



def write_output_csv_process(queue, out_file, total_num_entries):


    with tqdm(mininterval=0.5, unit="proteins", total=total_num_entries) as progress_bar:

        # TODO

        progress_bar.update() # update progress





if __name__ == "__main__":

    _, graph_gen_args, _ = parse_args()


    entry_queue = Queue(100)
    statistics_queue = Queue()








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
    graph_gen = [Process(target=generate_graph_consumer, args=(entry_queue, statistics_queue), kwargs=graph_gen_args) for _ in range(2)]
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
            entry = statistics_queue.get(timeout=3)
            # (prot, count) = prot_variation_queue.get()
            # t.append( entry)

        except Exception:

            pass




    # Nodes Count: 71031227
    # Edges Count: 77997502

    # Optimization is possible here!
    # We could concat Nodes together as long as there is only one  in and out degree
    # This optimization can happen before!!! doing the weighting (we can use the topological sort for this)

