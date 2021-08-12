import json

import cassandra
from cassandra.cluster import Cluster
from Bio.SwissProt import FeatureLocation, FeatureTable

from protgraph.export.abstract_exporter import AExporter

class Cassandra(AExporter):
    """
        TODO
    """

    def start_up(self, **kwargs):



        self.cluster = Cluster(["127.0.0.1"], port=9042)
        try:
            self.session = self.cluster.connect("graphs")
        except:
            try:
                self.session = self.cluster.connect()
                self.session.execute("CREATE KEYSPACE IF NOT EXISTS tttt WITH REPLICATION = { 'class' : 'SimpleStrategy', 'replication_factor' : 0 };")
                self.session.set_keyspace("keyspace")
            except Exception as e:
                raise Exception(e)  # TODO DL


    a = 1-3



    def export(self, prot_graph, _):
        # Export the protein
        self._export(prot_graph)

    def tear_down(self):
        # Close the connection to mysql
        try:
            self.cursor.close()  # Close cursor
            self.conn.close()  # Close connection
        except Exception as e:
            print("Connection to mysql could not be closed. (Reason: {})".format(str(e)))


    def _get_node_edge_attrs(self, node_edge_attrs, key_list):
        """ Get values of nodes/edges, returning None if not present """
        attrs_l = [None]*len(key_list)
        # Return list of node attrs
        for idx, ele in enumerate(key_list):
            if ele in node_edge_attrs:
                if ele == "qualifiers":
                    # Special Case for qualifiers here we do JSON!
                    attrs_l[idx] = json.dumps(self._get_attributes(node_edge_attrs[ele]))
                else:
                    attrs_l[idx] = node_edge_attrs[ele]
        return attrs_l

    def _get_attributes(self, attrs):
        """ Convert qualifiers objects into JSON-Serializable objects """
        if isinstance(attrs, list):
            return [self._get_attributes(x) for x in attrs]
        elif isinstance(attrs, dict):
            return {self._get_attributes(x): self._get_attributes(y) for x, y in attrs.items()}
        elif isinstance(attrs, FeatureLocation):
            return [attrs.nofuzzy_start, attrs.nofuzzy_end]
        elif isinstance(attrs, FeatureTable):
            return {self._get_attributes(x): self._get_attributes(y) for x, y in attrs.__dict__.items()}
        else:
            return attrs
