import redis
from Bio.SeqFeature import FeatureLocation
from Bio.SwissProt import FeatureTable
from redisgraph import Edge, Graph, Node

from protgraph.export.abstract_exporter import AExporter
from protgraph.graph_collapse_edges import Or


class RedisGraph(AExporter):
    """ A simple exporter to the module RedisGraph using RediGraph.py """

    def start_up(self, **kwargs):
        # Here we generate a connection to the redis server
        # which should have the module RedisGraph
        self.host = kwargs["redisgraph_host"]  # Host
        self.port = kwargs["redisgraph_port"]  # Port
        self.graph_name = kwargs["redisgraph_graph"]  # Name where we add graphs

        # Initialize connection
        try:
            self.conn = redis.Redis(host=self.host, port=self.port)
        except Exception as e:
            raise Exception("Could not establish a connection to redis", e)

    def export(self, prot_graph, _):
        # Export The protein
        self._export(prot_graph)

    def tear_down(self):
        # Close the connection to RedisGraph
        try:
            self.conn.close()
        except Exception as e:
            print("Connection to RedisGraph could not be closed. (Reason: {})".format(str(e)))

    def _export(self, prot_graph):
        """ Actual Export, simply generate the corresponding nodes and edges and commit it """
        # Create graph
        redis_graph = Graph(self.graph_name, self.conn)

        # Create nodes
        redis_nodes = [
            Node(label="node", properties=self._get_attributes_string(x.attributes()))
            for x in prot_graph.vs[:]
        ]
        for x in redis_nodes:
            redis_graph.add_node(x)

        # Create edges
        redis_edges = [
            Edge(
                redis_nodes[x.source], "edge", redis_nodes[x.target],
                properties=self._get_attributes_string(x.attributes())
            )
            for x in prot_graph.es[:]
        ]
        for x in redis_edges:
            redis_graph.add_edge(x)

        # Commit graph to the server
        try:
            redis_graph.commit()
        except Exception as e:
            raise Exception("Could not add Protein '{}' to RedisGraph! (Reason: {})"
                            .format(prot_graph.vs["accession"][0], str(e)))
            # TODO do we start again and retry once?

    def _get_attributes_string(self, attrs):
        """ We need to use this, since properties can only have depth == 1 (due to CYPHER?! / RedisGraph.py!?) """
        if attrs is None:
            return "NULL"
        elif isinstance(attrs, list):
            return "[" + ",".join([str(x) for x in attrs]) + "]"
        elif isinstance(attrs, dict):
            return {self._get_attributes_string(x): self._get_attributes_string(y) for x, y in attrs.items()}
        else:
            return attrs

    def _get_attributes(self, attrs):
        """ UNUSED: Alternative conversion of properties in more complex datatype. """
        if attrs is None:
            return "NULL"
        if isinstance(attrs, Or):
            return {"or": [self._get_attributes(x) for x in attrs]}
        elif isinstance(attrs, list):
            return [self._get_attributes(x) for x in attrs]
        elif isinstance(attrs, dict):
            return {self._get_attributes(x): self._get_attributes(y) for x, y in attrs.items()}
        elif isinstance(attrs, FeatureLocation):
            return [attrs.nofuzzy_start, attrs.nofuzzy_end]
        elif isinstance(attrs, FeatureTable):
            return {self._get_attributes(x): self._get_attributes(y) for x, y in attrs.__dict__.items()}
        else:
            return attrs
