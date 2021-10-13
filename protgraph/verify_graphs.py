def verify_graph(graph):

    _check_dag(graph)

    _check_direct_edge_to_end(graph)

    _check_parallel_edges(graph)

    # Check degrees!
    _check_degree(graph)
    _check_indegree(graph)
    _check_outdegree(graph)


def _check_direct_edge_to_end(graph):
    # Get start and end point
    [__start_node__] = graph.vs.select(aminoacid="__start__")
    [__stop_node__] = graph.vs.select(aminoacid="__end__")

    if __stop_node__.index in [x.index for x in __start_node__.neighbors(mode="OUT")]:
        print(
            "Protein {} has a direct edge from start to end! "
            "This is a case where a null protein/peptide is represented. This "
            "should not be valid!".format(graph.vs[0]["accession"])
        )


def _check_degree(graph):
    if len([x for x in graph.degree() if x == 0]) != 0:
        print(
            "Degree of Protein {} is for some vertices set at 0! "
            "There might be some not connected Subgraphs!".format(graph.vs[0]["accession"])
        )


def _check_indegree(graph):
    if len([x for x in graph.indegree() if x == 0]) != 1:
        print(
            "InDegree of Protein {} is not 1 (actual {})! "
            "There might be some not connected Subgraphs!".format(
                graph.vs[0]["accession"],
                len([x for x in graph.indegree() if x == 0])
            )
        )


def _check_outdegree(graph):
    if len([x for x in graph.outdegree() if x == 0]) != 1:
        print(
            "OutDegree of Protein {} is not 1 (actual {})! "
            "There might be some not connected Subgraphs!".format(
                graph.vs[0]["accession"],
                len([x for x in graph.outdegree() if x == 0])
            )
        )


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
