from contextlib import ContextDecorator

from protgraph.export.csv import CSV
from protgraph.export.dot import Dot
from protgraph.export.gml import GML
from protgraph.export.graphml import GraphML
from protgraph.export.gremlin import Gremlin
from protgraph.export.large_csv import Large_CSV
from protgraph.export.mysql import MySQL
from protgraph.export.peptides.pep_citus import PepCitus
from protgraph.export.peptides.pep_fasta import PepFasta
from protgraph.export.peptides.pep_mysql import PepMySQL
from protgraph.export.peptides.pep_postgres import PepPostgres
from protgraph.export.pickle import Pickle
from protgraph.export.postgres import Postgres
from protgraph.export.redisgraph import RedisGraph
from protgraph.export.cassandra import Cassandra


class Exporters(ContextDecorator):
    """
    A class containing a list of all exporters
    This class can be closed afterwards -> call it so that the exporters can tear down!
    """

    def __init__(self, **kwargs):
        """ Set a list of enabled exporters here. """
        self.export_classes = []

        if kwargs["export_dot"]:
            self.export_classes.append(Dot())
        if kwargs["export_csv"]:
            self.export_classes.append(CSV())
        if kwargs["export_large_csv"]:
            self.export_classes.append(Large_CSV())
        if kwargs["export_graphml"]:
            self.export_classes.append(GraphML())
        if kwargs["export_gml"]:
            self.export_classes.append(GML())
        if kwargs["export_pickle"]:
            self.export_classes.append(Pickle())
        if kwargs["export_redisgraph"]:
            self.export_classes.append(RedisGraph())
        if kwargs["export_postgres"]:
            self.export_classes.append(Postgres())
        if kwargs["export_mysql"]:
            self.export_classes.append(MySQL())
        if kwargs["export_gremlin"]:
            self.export_classes.append(Gremlin())
        if kwargs["export_peptide_postgres"]:
            self.export_classes.append(PepPostgres())
        if kwargs["export_peptide_mysql"]:
            self.export_classes.append(PepMySQL())
        if kwargs["export_peptide_fasta"]:
            self.export_classes.append(PepFasta())
        if kwargs["export_peptide_citus"]:
            self.export_classes.append(PepCitus())
        if kwargs["export_cassandra"] or True: # TODO remove
            self.export_classes.append(Cassandra())

        # Also start up all exporters
        for ec in self.export_classes:
            ec.start_up(**kwargs)

    def export_graph(self, prot_graph, out_queue):
        """ Mapping to export a protein graph to all exporters """
        for ec in self.export_classes:
            ec.export(prot_graph, out_queue)

    def close(self):
        """ Tear down all available exporters """
        for ec in self.export_classes:
            ec.tear_down()
