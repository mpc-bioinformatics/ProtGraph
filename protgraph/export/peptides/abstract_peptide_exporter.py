from abc import abstractmethod

import networkx

from protgraph.export.abstract_exporter import AExporter


class APeptideExporter(AExporter):
    """
    This is a abstract exporter for peptides in a graph. It can implement
    an export funtionality to a folder / database or more.
    """

    def __init__(self):
        """ Initialize needed parameters with None """
        self.skip_x = None  # Skip Peptides, which contain X?
        self.peptide_min_length = None  # Minimum peptide length
        self.max_miscleavages = None  # Maximum number of miscleavages in a peptide
        self.use_igraph = None  # Use Igraph? (or networkX?)
        self.peptide_max_length = None  # Maximum peptide length to limit possibilites. None to consider all.
        self.batch_size = None  # Batch size of peptides which will be processed at once. (list length)
        self.mass_factor = None  # Mass factor for only allowing min/max weight peptides
        self.min_weight = None  # The minimum weight of the peptide in Da
        self.max_weight = None  # The maximum weight of the peptide in Da

    def _set_up_taversal(
        self, skip_x, peptide_min_length, max_miscleavages, use_igraph,
        peptide_max_length, batch_size, min_weight, max_weight
    ):
        """
        Set parameters by dedicated function (preferably in start_up)
        This method is seperately since it may be needed that we want finer control in export at once
        E.G.: PGExport with 5-60 AAs, SQLiteExport with 5-50 etc...
        """
        self.skip_x = skip_x
        self.peptide_min_length = peptide_min_length
        self.max_miscleavages = max_miscleavages
        self.use_igraph = use_igraph
        self.peptide_max_length = peptide_max_length
        self.batch_size = batch_size
        self.min_weight = min_weight
        self.max_weight = max_weight

    def start_up(self, **kwargs):
        self._set_up_taversal(
            kwargs["pep_skip_x"],
            kwargs["pep_min_pep_length"],
            kwargs["pep_miscleavages"],
            kwargs["pep_use_igraph"],
            kwargs["pep_hops"],
            kwargs["pep_batch_size"],
            int(kwargs["pep_min_weight"] * kwargs["mass_dict_factor"]),
            int(kwargs["pep_max_weight"] * kwargs["mass_dict_factor"]),
        )

    @abstractmethod
    def export_peptides(self, prot_graph, l_path_nodes, l_path_edges, l_peptide, l_miscleavages, queue):
        """ Here goes the actual implementation of exporting a list of peptides into a file/folder/database etc. """
        pass

    def export(self, prot_graph, queue):
        """
        This abstract exporter implements the actual peptide traversal
        in a graph and calls the peptide exporter
        """
        # Get start and end node
        [__start_node__] = prot_graph.vs.select(aminoacid="__start__")
        [__stop_node__] = prot_graph.vs.select(aminoacid="__end__")

        # Set batch lists
        l_path = []
        l_edge_ids = []
        l_aas = []
        l_cleaved = []
        # Iterate over all peptides
        for path in self._get_peps(prot_graph, __start_node__, __stop_node__):

            # Get the actual peptide (concatenated aminoacids)
            aas = "".join(prot_graph.vs[path[1:-1]]["aminoacid"])

            # Get the edge ids from a path
            pairs = [(a, b) for a, b in zip(path, path[1:])]
            edge_ids = prot_graph.get_eids(pairs=pairs)

            # Skip Peptides, which contain an X
            if self.skip_x and "X" in aas:
                continue

            # Filter by peptide length
            if len(aas) < self.peptide_min_length:
                continue

            # Get number of cleaved edges
            if "cleaved" in prot_graph.es[edge_ids[0]].attributes():
                cleaved = sum(filter(None, prot_graph.es[edge_ids]["cleaved"]))
            else:
                cleaved = -1

            # And filter by miscleavages
            if self.max_miscleavages != -1:
                if cleaved > self.max_miscleavages:
                    continue

            if "mono_weight" in prot_graph.vs[0].attributes():
                if self.min_weight > 0:
                    w = sum(prot_graph.vs[path]["mono_weight"])
                    if self.min_weight > w:
                        continue
                if self.max_weight > 0:
                    w = sum(prot_graph.vs[path]["mono_weight"])
                    if w > self.max_weight:
                        continue

            # Append information to list
            l_path.append(path)
            l_edge_ids.append(edge_ids)
            l_aas.append(aas)
            l_cleaved.append(cleaved)

            if len(l_path) >= self.batch_size:
                # We export the list of peptides here and reset those lists afterwards
                self.export_peptides(prot_graph, l_path, l_edge_ids, l_aas, l_cleaved, queue)
                l_path = []
                l_edge_ids = []
                l_aas = []
                l_cleaved = []

        if len(l_path) > 0:
            # Special case, we might have some peptides left
            self.export_peptides(prot_graph, l_path, l_edge_ids, l_aas, l_cleaved, queue)

    def _get_peps(self, prot_graph, s, e):
        """ Get peptides depending on selected method """
        # OFFSET +1 since we have dedicated start and end nodes!
        # Except for -1, then we consider all paths
        if self.peptide_max_length < 0:
            if not self.use_igraph:
                cutoff = None
            else:
                cutoff = -1
        else:
            cutoff = self.peptide_max_length + 1

        if self.use_igraph:
            # This can consume lots of memory
            results = prot_graph.get_all_simple_paths(
                s.index,
                to=e.index,
                cutoff=cutoff
            )
            for r in results:
                yield r
        else:
            # This is a generator approach but is also considerably slower
            netx = prot_graph.to_networkx()
            yield from networkx.algorithms.simple_paths.all_simple_paths(
                netx,
                s.index,
                e.index,
                cutoff=cutoff
            )
