from contextlib import ContextDecorator

from export.dot import Dot
from export.graphml import GraphML
from export.gml import GML


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
        if kwargs["export_graphml"]:
            self.export_classes.append(GraphML())
        if kwargs["export_gml"]:
            self.export_classes.append(GML())

        # Also start up all exporters
        for ec in self.export_classes:
            ec.start_up(**kwargs)

    def export_graph(self, prot_graph):
        """ Mapping to export a protein graph to all exporters """
        for ec in self.export_classes:
            ec.export(prot_graph)

    def close(self):
        """ Tear down all available exporters """
        for ec in self.export_classes:
            ec.tear_down()
