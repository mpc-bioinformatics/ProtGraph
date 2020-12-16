from multiprocessing import Queue, Process


from read_embl import read_embl
from graph_generator import generate_graph
from variation_generator import get_next_variant
from database_writer import insert_to_database
import time


import psycopg2

if __name__ == "__main__":
    entry_queue = Queue(15000)
    graph_queue = Queue(15000)
    prot_variation_queue = Queue(15000)
    output_queue = Queue(15000)







    input_file = ["/home/luxdo/git/uniprot-filtered-organism Homo+sapiens+(Human)+[9606] .txt"] # TODO parameterize
    input_file_totals = [194237]

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
    graph_gen = [Process(target=generate_graph, args=(entry_queue, graph_queue,)) for _ in range(4)]
    prot_var_gen = [Process(target=get_next_variant, args=(graph_queue, prot_variation_queue,)) for _ in range(44)]
    database_writer = [Process(target=insert_to_database, args=(prot_variation_queue, output_queue,)) for _ in range(16)]
    


    # Start Processes in reverse!
    for p in database_writer:
        p.start()
    for p in prot_var_gen:
        p.start()
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
            entry = output_queue.get(timeout=60)
            # (prot, count) = prot_variation_queue.get()
            # t.append( entry)

        except Exception:

            pass




    # Nodes Count: 71031227
    # Edges Count: 77997502

    # Optimization is possible here!
    # We could concat Nodes together as long as there is only one  in and out degree
    # This optimization can happen before!!! doing the weighting (we can use the topological sort for this)
