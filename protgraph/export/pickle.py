from protgraph.export.generic_file_exporter import GenericFileExporter


class Pickle(GenericFileExporter):
    """ A simple GraphML exporter """

    def __init__(self):
        super(Pickle, self).__init__(
            lambda pg, path: pg.write_pickle(path + ".pickle")
        )
