import json

# import cassandra
from Bio.SeqFeature import FeatureLocation
from Bio.SwissProt import FeatureTable
# from cassandra.query import BatchStatement
from cassandra import InvalidRequest
from cassandra.cluster import Cluster

from protgraph.export.abstract_exporter import AExporter
# import cassandra.concurrent
from protgraph.graph_collapse_edges import Or


class Cassandra(AExporter):
    """
    This is a cassandra exporter, which is able to export nodes and edges
    into a cassandra cluster.
    """

    def start_up(self, **kwargs):
        # Set up host, port and the keyspace.
        host = kwargs["cassandra_host"]
        port = kwargs["cassandra_port"]
        keyspace = kwargs["cassandra_keyspace"]
        # the chunk size in cassandra is limited by kb, but we are referencing
        # number of elements with this parameter. So keep it sufficiently small.
        self.chunk_size = kwargs["cassandra_chunk_size"]

        # Initialize COnnection
        not_init = True
        while not_init:
            try:
                self.cluster = Cluster([host], port=port)
                self.session = self.cluster.connect(keyspace)
                not_init = False
            except Exception:
                try:
                    self.cluster = Cluster([host], port=port)
                    self.session = self.cluster.connect()
                    self.session.execute(
                        "CREATE KEYSPACE IF NOT EXISTS " + keyspace +
                        " WITH REPLICATION = { 'class' : 'SimpleStrategy', 'replication_factor' : 1 };"
                    )
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
        # mono_weight, avrg_weight, mono_weight_to_end, avrg_weight_to_end
        # Edges:
        # qualifiers (List),  cleaved
        try:
            # create nodes
            self.session.execute("""
                CREATE TABLE IF NOT EXISTS nodes (
                    id bigint PRIMARY KEY,
                    accession ascii,
                    aminoacid ascii,
                    position int,
                    isoform_accession ascii,
                    isoform_position int,
                    mono_weight {0},
                    mono_weight_to_end {0},
                    avrg_weight {0},
                    avrg_weight_to_end {0}
                );""".format("bigint" if kwargs["mass_dict_type"] is int else "decimal"))

        except Exception as e:
            print("Warning: Failed creating table 'nodes' (Reason: {})".format(str(e)))
        finally:
            self.nodes_keys = [
                "accession",
                "aminoacid",
                "position",
                "isoform_accession",
                "isoform_position",
                "mono_weight",
                "mono_weight_to_end",
                "avrg_weight",
                "avrg_weight_to_end",
            ]

        try:
            # Create edges
            self.session.execute("""
                create table if not exists edges (
                    id bigint PRIMARY KEY,
                    source bigint,
                    target bigint,
                    cleaved boolean,
                    qualifiers text
                );""")

        except Exception as e:
            print("Warning: Failed creating table 'edges' (Reason: {})".format(str(e)))
        finally:
            self.edges_keys = [
                "cleaved",
                "qualifiers"
            ]

        # Set Input of nodes
        self.prep_stmt_nodes = self.session.prepare(
            "INSERT INTO nodes (id," + ",".join(self.nodes_keys) + ") "
            "VALUES (?," + ",".join(["?"]*len(self.nodes_keys)) + ")"
        )
        self.prep_insert_nodes_lambda = lambda x: self.session.execute_async(self.prep_stmt_nodes, x)

        # Set Input of edges
        self.prep_stmt_edges = self.session.prepare(
            "INSERT INTO edges (id,source,target," + ",".join(self.edges_keys) + ") "
            "VALUES (?,?,?," + ",".join(["?"]*len(self.edges_keys)) + ")"
        )
        self.prep_insert_edges_lambda = lambda x: self.session.execute_async(self.prep_stmt_edges, x)

        self.nodes_counter = self.unique_id_gen(**kwargs)
        self.edges_counter = self.unique_id_gen(**kwargs)

    def export(self, prot_graph, _):
        # Export the protein
        # Here two versions exist:
        # 1. Prepared Statement (currently uncommented)
        # 2. Concurrent Example (may be slow)

        try:
            # Concurrent example of adding nodeas and edges:
            # prep_ins = [
            #     (
            #         self.prep_stmt_nodes,
            #         [next(self.nodes_counter), *self._get_node_edge_attrs(x.attributes(), self.nodes_keys)]
            #     )
            #     for x in prot_graph.vs[:]
            # ]
            # cassandra.concurrent.execute_concurrent(self.session, prep_ins)
            # prep_ins_edges = [
            #     (self.prep_stmt_edges, [
            #         next(self.edges_counter),
            #         prep_ins[x.source][1][0],
            #         prep_ins[x.target][1][0],
            #         *self._get_node_edge_attrs(x.attributes(), self.edges_keys)
            #     ])
            #     for x in prot_graph.es[:]
            # ]
            # cassandra.concurrent.execute_concurrent(self.session, prep_ins_edges)

            # Prepared Statement Example:
            # Add nodes
            nodes = [
                [next(self.nodes_counter), *self._get_node_edge_attrs(x.attributes(), self.nodes_keys)]
                for x in prot_graph.vs[:]
            ]
            # batch_nodes = BatchStatement()
            for n_chunk in self.chunked_iterable(nodes, self.chunk_size):  # Longest possible batch
                for n in n_chunk:
                    # batch_nodes.add(self.prep_stmt_nodes, n)
                    self.prep_insert_nodes_lambda(n)
                # self.session.execute(batch_nodes)
                # batch_nodes.clear()

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
            # batch_edges = BatchStatement()
            for e_chunk in self.chunked_iterable(edges, self.chunk_size):  # Longest possible batch
                for e in e_chunk:
                    # batch_nodes.add(self.prep_stmt_edges, e)
                    self.prep_insert_edges_lambda(e)
                # self.session.execute(batch_edges)
                # batch_edges.clear()

        except Exception as e:
            # TODO custom exception, this should not happen!
            # track information like key, which node/edge especially
            # retrieve the accession explicitly here too for debugging!
            print("ERROR adding Node/Edges of Protein {} into Cassandra. Please report this issue: {}".format(
                prot_graph.vs[0]["accession"], e
            ))
        except InvalidRequest as ir:
            print("Error Chunk size may be too large. Reason: {}".format(str(ir)))

    def tear_down(self):
        # Close the connection to mysql
        try:
            self.session.shutdown()
            self.cluster.shutdown()
        except Exception as e:
            print("Connection to Cassandra could not be shutdown. (Reason: {})".format(str(e)))

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
        if isinstance(attrs, Or):
            return {"or": [self._get_attributes(x) for x in attrs]}
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
