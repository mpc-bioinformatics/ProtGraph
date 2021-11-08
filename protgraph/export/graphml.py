from protgraph.export.generic_file_exporter import GenericFileExporter


class GraphML(GenericFileExporter):
    """ A simple GraphML exporter """

    def __init__(self):
        super(GraphML, self).__init__(
            lambda pg, path: pg.write_graphml(path + ".graphml")
        )
