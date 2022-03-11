import itertools
from abc import ABC, abstractmethod


class AExporter(ABC):
    """
    This is a abstract exporter for a protein graph. It can implement
    a export funtionality to a folder / database or more.
    """

    @abstractmethod
    def start_up(self, **kwargs):
        """
        This should implement the start up routine for a exporter.
        E.G. databases: Here a connection should be established.
        E.G. file Exporters: Here should the output folder be generated

        NOTE this is executed per Process! So make sure that this code can
        be called multiple times!
        """
        pass

    @abstractmethod
    def export(self, prot_graph, out_queue):
        """ Here goes the actual implementation for exporting! """
        pass

    @abstractmethod
    def tear_down(self):
        """
        Similar to the start_up routine, e.g. connections to databases should be closed here
        and opened file descriptors closed.!

        NOTE: This tear_down routine is also executed per Process. Make sure you do not close something,
        causing exceptions for other processes!
        """
        pass

    def unique_id_gen(self, **kwargs):
        val = kwargs["proc_id"]
        jump = kwargs["num_of_processes"]

        while True:
            yield val
            val += jump

    def chunked_iterable(self, iterable, size):
        """ Chunk down an iterable to a specific size """
        it = iter(iterable)
        while True:
            chunk = tuple(itertools.islice(it, size))
            if not chunk:
                break
            yield chunk
