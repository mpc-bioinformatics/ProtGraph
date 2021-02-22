def verify_graph(graph):

    _check_dag(graph)

    _check_parallel_edges(graph)


def _check_dag(graph):
    if not graph.is_dag:
        print("Protein {} is not a DAG!".format(graph.vs[0]["accession"]))


def _check_parallel_edges(graph):
    for x in range(0, graph.vcount()):
        k = graph.neighbors(x)
        k_set = set(k)
        if len(k) != len(k_set):
            print(
                "Protein '{}' has {} parallel edges (on Node {})".format(
                    graph.vs[0]["accession"],
                    len(k) - len(k_set),
                    x
                )
            )
