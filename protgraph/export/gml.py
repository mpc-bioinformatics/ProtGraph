from protgraph.export.generic_file_exporter import GenericFileExporter


class GML(GenericFileExporter):
    """ A simple GML (Graph Markup Language) exporter """

    def __init__(self):
        super(GML, self).__init__(
            lambda pg, path: pg.write_gml(path + ".gml")
        )
