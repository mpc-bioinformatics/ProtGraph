from abc import ABC, abstractmethod

class AExporter(ABC):
    ''' 
        This is a abstract exporter for a protein graph. It can implement
        a export funtionality to a folder / database or more. 
    '''

    @abstractmethod
    def start_up(self, **kwargs):
        '''
            This should implement the start up routine for a exporter.
            E.G. databases: Here a connection should be established.
            E.G. file Exporters: Here should the output folder be generated
            
            NOTE this is executed per Process! So make sure that this code can
            be called multiple times!
        '''
        pass

    @abstractmethod
    def export(self, prot_graph):
        ''' Here goes the actual implementation for exporting! '''
        pass

    @abstractmethod
    def tear_down(self):
        '''
            Similar to the start_up routine, e.g. connections to databases should be closed here
            and opened file descriptors closed.!

            NOTE: This tear_down routine is also executed per Process. Make sure you do not close something, 
            causing exceptions for other processes!
        '''
        pass
