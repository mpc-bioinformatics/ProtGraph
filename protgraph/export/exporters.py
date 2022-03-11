import importlib
from contextlib import ContextDecorator


class Exporters(ContextDecorator):
    """
    A class containing a list of all exporters
    This class can be closed afterwards -> call it so that the exporters can tear down!
    """

    def __enter__(self):
        """ Enter for ContextManager """
        return self

    def __init__(self, **kwargs):
        """ Set a list of enabled exporters here. """
        self.export_classes = []

        # Go through all implemented exporters and add if needed
        for entries in self.__available_exporters(**kwargs):
            self.__add_generic(self.export_classes, *entries)

        # Also start up all exporters
        for ec in self.export_classes:
            ec.start_up(**kwargs)

    def export_graph(self, prot_graph, out_queue):
        """ Mapping to export a protein graph to all exporters """
        for ec in self.export_classes:
            ec.export(prot_graph, out_queue)

    def __exit__(self, *exc):
        """ Tear down all available exporters """
        for ec in self.export_classes:
            ec.tear_down()

    # Actual supported Exporters are below
    def __add_generic(self, exporters, flag, export_from, export_class):
        """ Load exporting class dynamically """
        if flag:
            try:
                module = importlib.import_module(export_from)  # from x.y.z
                clazz = getattr(module, export_class)  # import class
                exporters.append(clazz())
            except ImportError as e:
                print("Skipping export into {export_name}! This export misses the dependency '{dep}'".format(
                    export_name=export_class,
                    dep=e.name
                ))

    def __available_exporters(self, **kwargs):
        """ List of available exporters """
        return [  # (FLAG, MODULE_LOC, CLASS_NAME)
            # Exporter to single graph files
            (kwargs["export_dot"], "protgraph.export.dot", "Dot"),
            (kwargs["export_csv"], "protgraph.export.csv", "CSV"),
            (kwargs["export_graphml"], "protgraph.export.graphml", "GraphML"),
            (kwargs["export_gml"], "protgraph.export.gml", "GML"),
            (kwargs["export_pickle"], "protgraph.export.pickle", "Pickle"),
            (kwargs["export_pcsr"], "protgraph.export.pcsr", "PCSR"),
            (kwargs["export_binary_pcsr"], "protgraph.export.binary_pcsr", "BinaryPCSR"),

            # Exporter to ONE large file
            (kwargs["export_large_csv"], "protgraph.export.large_csv", "LargeCSV"),
            (kwargs["export_large_pcsr"], "protgraph.export.large_pcsr", "LargePCSR"),
            (kwargs["export_large_binary_pcsr"], "protgraph.export.large_binary_pcsr", "LargeBinaryPCSR"),

            # Exporter to "databases" (setup required)
            (kwargs["export_postgres"], "protgraph.export.postgres", "Postgres"),
            (kwargs["export_redisgraph"], "protgraph.export.redisgraph", "RedisGraph"),
            (kwargs["export_mysql"], "protgraph.export.mysql", "MySQL"),
            (kwargs["export_gremlin"], "protgraph.export.gremlin", "Gremlin"),
            (kwargs["export_cassandra"], "protgraph.export.cassandra", "Cassandra"),

            # Peptide Exporter to "databases" (setup required)
            (kwargs["export_peptide_citus"], "protgraph.export.peptides.pep_citus", "PepCitus"),
            (kwargs["export_peptide_postgres"], "protgraph.export.peptides.pep_postgres", "PepPostgres"),
            (kwargs["export_peptide_mysql"], "protgraph.export.peptides.pep_mysql", "PepMySQL"),

            # Peptide Exporter to local filesystem (no setup required)
            (kwargs["export_peptide_fasta"], "protgraph.export.peptides.pep_fasta", "PepFasta"),
            (kwargs["export_peptide_trie"], "protgraph.export.peptides.pep_trie", "PepTrie"),
            (kwargs["export_peptide_sqlite"], "protgraph.export.peptides.pep_sqlite", "PepSQLite"),
        ]
