from multiprocessing import Queue, Process


from read_embl import read_embl
from graph_generator import generate_graph
from variation_generator import get_next_variant


if __name__ == "__main__":
    entry_queue = Queue(500)
    graph_queue = Queue(500)
    prot_variation_queue = Queue(500)


    input_file = "/home/luxii/git/PROVAR/uniprot-filtered-organism Homo+sapiens+(Human)+[9606] .txt" # TODO parameterize

    # Create Processes
    entry_reader = Process(target=read_embl, args=(input_file, entry_queue,))
    graph_gen = [Process(target=generate_graph, args=(entry_queue, graph_queue,)) for _ in range(1)]
    prot_var_gen = [Process(target=get_next_variant, args=(graph_queue, prot_variation_queue,)) for _ in range(3)]
    


    # Start Processes in reverse!
    for p in prot_var_gen:
        p.start()
    for p in graph_gen:
        p.start()
    entry_reader.start()






    t = []
    while True: 
        
        try:
            (prot, count) = prot_variation_queue.get(timeout=1)
            t.append( (prot, count) )

        except Exception:
            pass
