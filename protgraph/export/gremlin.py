import random
from time import sleep

from gremlin_python.driver.driver_remote_connection import \
    DriverRemoteConnection
from gremlin_python.process.anonymous_traversal import traversal
from gremlin_python.statics import FloatType, LongType

from protgraph.export.abstract_exporter import AExporter


class Gremlin(AExporter):
    """
    A Gremlin exporter (to various databases)

    NOTE: This exporter is not finisehd and may not be further developed.
    """

    def start_up(self, **kwargs):
        # it can cause exceptions when starting all at once!
        # This solves it currently... ( TODO why does it not allow multiple connections consistently?)
        sleep(random.randint(1, 15))

        # This generates a connection to the gremlin server
        self.conn = DriverRemoteConnection(kwargs["gremlin_url"], kwargs["gremlin_traversal_source"], pool_size=1)
        self.remote_graph = traversal().withRemote(self.conn)

    def export(self, prot_graph, _):
        try:
            # Non Bulk approach
            nodes = [
                self._set_properties(self.remote_graph.addV(), list(x.attributes().items())).next()
                for x in prot_graph.vs[:]
            ]
            _ = [
                self._set_properties(
                    self.remote_graph.addE("e").from_(nodes[x.source]).to(nodes[x.target]),
                    list(x.attributes().items())
                ).next()
                for x in prot_graph.es[:]
            ]

        except Exception as e:
            raise e

        # Bulk approach, which works only partially (due to high processing time on server side!)
        # try:
        #     rg = self.remote_graph
        #     for k in prot_graph.vs:
        #         rg = self._set_properties(rg.addV(), list(k.attributes().items())).as_(str(k.index))
        #     for k in prot_graph.es:
        #         rg = self._set_properties(
        #             rg.addE("e").from_(str(k.source)).to(str(k.target)),
        #             list(k.attributes().items())
        #         )
        #     rg.next()
        # except Exception as e:
        #     raise e

    def _set_properties(self, rg, prot_graph_attrs):
        """ set the properties of nodes or edges recursively """
        if len(prot_graph_attrs) == 0:
            return rg
        key, value = prot_graph_attrs.pop()
        if key == "position" and value is None:
            return self._set_properties(rg, prot_graph_attrs).property(key, -1)
        if (key == "isoform_position" or key == "isoform_accession" or key == "cleaved" or key == "qualifiers") \
           and value is None:
            return self._set_properties(rg, prot_graph_attrs).property(key, "")
        if key == "qualifiers":
            # Special case, do a json dump to string
            return self._set_properties(rg, prot_graph_attrs).property(key, "")
        if (key == "avrg_weight" or key == "mono_weight" or key == "avrg_weight_to_end" or key == "mono_weight_to_end"):
            # Special case for weights
            if type(value) is int:
                return self._set_properties(rg, prot_graph_attrs).property(key, LongType(value))
            return self._set_properties(rg, prot_graph_attrs).property(key, FloatType(value))

        return self._set_properties(rg, prot_graph_attrs).property(key, value)

    def tear_down(self):
        # Close the connection to the gremlin server
        self.conn.close()
