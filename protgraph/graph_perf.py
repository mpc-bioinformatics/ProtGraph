import timeit

import igraph
import networkit as nk
import networkx

igraph_graph = igraph.load("exported_graphs/P/0/4/6/3/7.pickle")
# igraph_graph = igraph.load("exported_graphs/P/2/0/7/2/9.pickle")

CUTOFF= 4
[__start_node__] = igraph_graph.vs.select(aminoacid="__start__")
[__stop_node__] = igraph_graph.vs.select(aminoacid="__end__")
start = __start_node__.index
end = __stop_node__.index


### Method 3
def _find_all_paths(graph, start, end, cutoff):  # IGRAPH IMPL
    # Maybe in C or CPP ? Optimize this, it is not fast enough!
    # Generator for all possible paths
    path  = []
    queue = [(start, end, path)]
    while queue:
        start, end, path = queue.pop()
        if len(path) > cutoff: # TODO remove this MODIFICATION of the all path algorithm here!
            continue
        path = path + [start]

        for node in set(graph.neighbors(start,mode='OUT')).difference(path):
            queue.append((node, end, path))

        if start == end and len(path) > 0:              
            yield path

### method 4
def _find_all_paths_networkx(graph, start, end, cutoff):  # IGRAPH IMPL
    # Maybe in C or CPP ? Optimize this, it is not fast enough!
    # Generator for all possible paths
    path  = []
    queue = [(start, end, path)]
    while queue:
        start, end, path = queue.pop()
        if len(path) > cutoff: # TODO remove this MODIFICATION of the all path algorithm here!
            continue
        path = path + [start]

        for node in set(graph.successors(start)).difference(path):
            queue.append((node, end, path))

        if start == end and len(path) > 0:              
            yield path


def method_1():
    results = igraph_graph.get_all_simple_paths(__start_node__, to=__stop_node__, cutoff=CUTOFF)
    k = 0
    for x in results:
        k += 1


def method_2():
    netx = igraph_graph.to_networkx()
    k = 0
    for x in networkx.algorithms.simple_paths.all_simple_paths(netx, start, end, cutoff=CUTOFF):
        # RETURNS REPEATING RESULTS!
        k += 1


def method_3():
    k = 0
    for x in _find_all_paths(igraph_graph, __start_node__.index, __stop_node__.index, cutoff=CUTOFF):
        k += 1

def method_4():
    netx = igraph_graph.to_networkx()
    k = 0
    for x in _find_all_paths_networkx(igraph_graph, start, end, cutoff=CUTOFF):
        k += 1


def method_5():
    netx = igraph_graph.to_networkx()
    netkit = nk.nxadapter.nx2nk(netx)
    t = nk.distance.AllSimplePaths(netkit, start, end, cutoff=CUTOFF)
    t.run()
    k = 0
    def temp(j):
        k + 1
    t.forAllSimplePaths(temp)



print("igraph get_all_simple_paths", timeit.repeat(stmt=method_1, number=100))
print("networkx all_simple_paths  ", timeit.repeat(stmt=method_2, number=100))
print("igraph manual_traversal    ", timeit.repeat(stmt=method_3, number=100))
print("networkx manual_traversal  ", timeit.repeat(stmt=method_4, number=100))
print("networkit AllSimplePaths   ", timeit.repeat(stmt=method_5, number=100))