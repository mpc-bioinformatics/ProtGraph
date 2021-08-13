import json
import itertools

from cassandra.cluster import Cluster
from cassandra.query import BatchStatement
from Bio.SwissProt import FeatureLocation, FeatureTable

from protgraph.export.abstract_exporter import AExporter
from cassandra import InvalidRequest

class Cassandra(AExporter):
    """
        TODO
    """

    def start_up(self, **kwargs):

        host = kwargs["cassandra_host"]
        port = kwargs["cassandra_port"]
        keyspace = kwargs["cassandra_keyspace"]
        self.chunk_size = kwargs["cassandra_chunk_size"]

        not_init = True
        while not_init:
            try:
                self.cluster = Cluster([host], port=port)
                self.session = self.cluster.connect(keyspace)
                not_init = False
            except:
                try:
                    self.cluster = Cluster([host], port=port)
                    self.session = self.cluster.connect()
                    self.session.execute("CREATE KEYSPACE IF NOT EXISTS " + keyspace + " WITH REPLICATION = { 'class' : 'SimpleStrategy', 'replication_factor' : 1 };")
                    self.session.set_keyspace(keyspace)
                    not_init = False
                except Exception as e:
                    print("Warning failed to create or change to keyspace. Reason: {}".format(str(e)))

        # Create Nodes and Edges tables in Cassandra
        # also generates the prepared statements
        self._create_tables(**kwargs)


    def _create_tables(self, **kwargs):
        """ Create the nodes and edges tables """
        # All currently used keys:
        # Nodes:
        # accession, aminoacid, position, isoform_accession, isoform_position
        # Edges:
        # qualifiers (List), mono_weight, avrg_weight, mono_weight_to_end, avrg_weight_to_end, cleaved
        try:
            # create nodes
            self.session.execute("""
                CREATE TABLE IF NOT EXISTS nodes (
                    id bigint PRIMARY KEY,
                    accession ascii,
                    aminoacid ascii,
                    position int,
                    isoform_accession ascii,
                    isoform_position int
                );""")

        except Exception as e:
            print("Warning: Failed creating table 'nodes' (Reason: {})".format(str(e)))
        finally:
            self.nodes_keys = [
                "accession",
                "aminoacid",
                "position",
                "isoform_accession",
                "isoform_position"
            ]

        try:
            # Create edges
            self.session.execute("""
                create table if not exists edges (
                    id bigint PRIMARY KEY,
                    source bigint,
                    target bigint,
                    cleaved boolean,
                    mono_weight {0},
                    mono_weight_to_end {0},
                    avrg_weight {0},
                    avrg_weight_to_end {0},
                    qualifiers text
                );""".format("bigint" if kwargs["mass_dict_type"] is int else "decimal"))

        except Exception as e:
            print("Warning: Failed creating table 'edges' (Reason: {})".format(str(e)))
        finally:
            self.edges_keys = [
                "cleaved",
                "mono_weight",
                "mono_weight_to_end",
                "avrg_weight",
                "avrg_weight_to_end",
                "qualifiers"
            ]

        # Set Input of nodes
        self.prep_stmt_nodes = self.session.prepare(
            "INSERT INTO nodes (id," + ",".join(self.nodes_keys) + ") VALUES (?," + ",".join(["?"]*len(self.nodes_keys)) + ")"
        )
        self.prep_insert_nodes_lambda = lambda x: self.session.execute(self.prep_stmt_nodes, x)


        # Set Input of edges
        self.prep_stmt_edges = self.session.prepare(
            "INSERT INTO edges (id,source,target," + ",".join(self.edges_keys) + ") VALUES (?,?,?," + ",".join(["?"]*len(self.edges_keys)) + ")"
        )
        self.prep_insert_edges_lambda = lambda x: self.session.execute(self.prep_stmt_edges, x)

        self.nodes_counter = self.unique_id_gen(**kwargs)
        self.edges_counter = self.unique_id_gen(**kwargs)
        


    def export(self, prot_graph, _):
        # Export the protein
        # Add nodes
        try:
            nodes = [[next(self.nodes_counter), *self._get_node_edge_attrs(x.attributes(), self.nodes_keys)] for x in prot_graph.vs[:]]
            for n_chunk in self._chunked_iterable(nodes, self.chunk_size):  # Longest possible batch
                batch_nodes = BatchStatement()
                for n in n_chunk:
                    batch_nodes.add(self.prep_stmt_nodes, n)
                    # self.prep_insert_nodes_lambda(n)
                self.session.execute(batch_nodes)

            # Add edges
            edges = [
                [
                    next(self.edges_counter), 
                    nodes[x.source][0],
                    nodes[x.target][0],
                    *self._get_node_edge_attrs(x.attributes(), self.edges_keys)
                ] 
                for x in prot_graph.es[:]
            ]
            for e_chunk in self._chunked_iterable(edges, self.chunk_size):  # Longest possible batch
                batch_edges = BatchStatement()
                for e in e_chunk:
                    batch_nodes.add(self.prep_stmt_edges, e)
                    # self.prep_insert_edges_lambda(e)
                self.session.execute(batch_edges)


        except Exception as e:
            print("help!")
        except InvalidRequest as ir:
            print("Error Chunk size may be too large. Reason: {}".format(str(ir)))


    def tear_down(self):
        # Close the connection to mysql
        try:
            self.session.shurdown()
            self.cluster.shutdown()
        except Exception as e:
            print("Connection to Cassandra could not be shutdown. (Reason: {})".format(str(e)))


    def _chunked_iterable(self, iterable, size):
        it = iter(iterable)
        while True:
            chunk = tuple(itertools.islice(it, size))
            if not chunk:
                break
            yield chunk


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
