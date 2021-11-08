from protgraph.export.generic_file_exporter import GenericFileExporter


class Dot(GenericFileExporter):
    """ A simple Dot exporter """

    def __init__(self):
        super(Dot, self).__init__(
            lambda pg, path: pg.write_dot(path + ".dot")
        )
