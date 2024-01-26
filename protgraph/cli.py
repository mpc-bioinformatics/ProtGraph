import os
import string
from argparse import ArgumentTypeError
from itertools import product


def check_if_file_exists(s: str):
    """ checks if a file exists. If not: raise Exception """
    if os.path.isfile(s):
        return s
    else:
        raise Exception("File '{}' does not exists".format(s))


def add_main_args(parser):
    # Needed Arguments for parsing (and other general/global arguments)
    parser.add_argument(
        "files", type=check_if_file_exists, nargs="+",
        help="Files containing the Swissprot/EMBL-Entries (either in .dat or .txt)"
    )
    parser.add_argument(
        "--num_of_entries", "-n", type=int, default=None,
        help="Number of entries across all files (summed). if given, it will an estimation of the runtime"
    )
    parser.add_argument(
        "--exclude_accessions", "-exclude", type=str, default=None,
        help="A csv file only containing accessions in the first row which should be excluded for processsing."
        " Setting this value may reduce the reading performance and therefore the throughput performance overall."
    )
    parser.add_argument(
        "--num_of_processes", "-np", type=int, default=None,
        help="The number of processes used to process entries. Each process can process an entry individually. "
        "The default value is 'cores - 1', since one process is reserved for reading entries"
    )
    parser.add_argument(
        "--output_csv", "-o", default=os.path.join(os.getcwd(), "protein_graph_statistics.csv"),
        type=str,
        help="Set the output file, which will contain information about the protein graphs (in csv) NOTE: "
        "It will write to 'protein_graph_statistics.csv' and overwrite if such a file exists. Default is "
        "set to the current working directory"
    )

    # Available headers
    avail_headers = [
        # Basic Information
        "accession", "gene_id",
        # Num of Fts
        "num_var_seq", "num_init_met", "num_signal", "num_variant", "num_mutagen",
        "num_conflict", "num_peptide", "num_propep", "num_chain", "num_nodes", "num_edges",
        # Statistics
        "num_paths", "list_paths_miscleavages", "list_paths_hops",
        "list_paths_variant", "list_paths_mutagen", "list_paths_conflict",
        # Protein description (at the end, since it can be very lengthy)
        "protein_description",
        ]
    # Add list functionality to CLI

    def _check_if_in_list(input: str):
        """
        Check if entries, seperated by ',' are all in the list
        Returns the list as by user
        """
        header_l = []
        for i in input.split(","):
            if i not in avail_headers:
                raise ArgumentTypeError(
                    "The header-Entry '{}' does not exists.".format(i)
                )
            header_l.append(i)
        return header_l
    parser.add_argument(
        "--output_csv_layout", type=_check_if_in_list, action="store", default=avail_headers,
        help="Set the csv layout of the generated graph statistics. You can choose from: " + ",".join(avail_headers) +
        " (this is also the default order.) You can specify your own order or exclude headers,"
        " by using a comma-seperated list."
    )

    parser.add_argument(
        "--no_description", "-no_desc", default=False, action="store_true",
        help="Set this flag to not include the protein descriptions into the output statistics file. "
        "This reduces the size of the output statistics if set."
    )


def add_graph_generation(group):
    # Available features
    avail_fts = [
        "ALL", "NONE", "INIT_MET", "VARIANT", "VAR_SEQ", "SIGNAL", "MUTAGEN", "CONFLICT", "PEPTIDE", "PROPEP", "CHAIN"
    ]
    group.add_argument(
        "--feature_table", "-ft", choices=avail_fts, type=str.upper, action="append",
        help="Set the features which should be added to the graph. Those feature_tables correspond to"
        " the actual FT name provided in SwissProt-EMBL. Default is set to use all available features"
        " ProtGraph can parse. Currently parsable features are: " + ", ".join(avail_fts) +
        ". Use it as follows to only select specific ones: '-ft INIT_MET -ft SIGNAL' or '-ft NONE' to use none."
    )

    # Add replace funcitonality to CLI
    def _replace_syntax(input: str):
        """
        Check if the replacement syntax is in form.
        Returns a tuple (s, t) where the amino acids gets replaced from s -> t[0], t[1], ...
        (similar to the MUTAGEN Syntax in SP-EMBL!)
        """
        if "->" not in input:
            raise ArgumentTypeError(
                "'->' is not set in replacement rule! Did you quoted the replacement as here: \"X->Y\"?"
            )
        s, t = input.split("->", 1)

        s = s.strip().upper()
        if not s.isalpha() or len(s) != 1:
            raise ArgumentTypeError(
                "The amino acid which gets replaced can only be set to: [A-Z] (1 letter)! Found: '{}'".format(s)
            )

        t = [x.strip().upper() for x in t.split(",")]
        for x in t:
            if not x.isalpha() or len(s) != 1:
                raise ArgumentTypeError(
                    "The amino acids to replace to can only be set to: [A-Z] (1 letter)! Found '{}'".format(x)
                )

        return s, t
    group.add_argument(
        "--replace_aa", "-raa", type=_replace_syntax, action="append",
        help="Substitute amino acids in graphs by other amino acids. This could be useful to replace"
        " e.g. 'J' with 'I' and 'L'. This parameter can be then provided as: 'J->I,L'. Multiple replacements"
        " are allowed and are executed one after another. NOTE: only from ONE amino acid multiple amino acids can"
        " be substituted. So only the format: 'A->B[,C]*' is allowed!"
    )

    def _parse_mod(input: str):
        """
        Parse Modification
        """
        if ":" not in input:
            raise ArgumentTypeError(
                "':' is not set in modification rule! Did you quoted the replacement as here: \"M:15.994915\"?"
            )
        s, t = input.split(":", 1)

        s = s.strip().upper()
        if s not in [
                    "CPEPTERM", "NPEPTERM", "NPROTERM", "CPROTERM",
                    *[x+y for x, y in list(product(["NPEP", "CPEP", ""], string.ascii_uppercase))]
                ]:  # Checks for C-/N-Term for Peptide/Protein as well as for [A-Z], NPEP[A-Z] and CPEP[A-Z]
            raise ArgumentTypeError(
                "The amino acid which gets replaced can only be set to: [A-Z] (1 letter)! Found: '{}'".format(s)
            )

        try:
            t = float(t)
        except Exception():
            raise ArgumentTypeError(
                "The DeltaMass is not a numeric number. Received: '{}'".format(t)
            )

        return s, t
    group.add_argument(
        "--fixed_mod", "-fm", type=_parse_mod, action="append",
        help="Apply a fixed modification on a special aminoacid."
        " You can apply multiple fix modifications BUT only one modification per aminoacid is currently allowed."
        " The form should be '<AminoAcid>:<DeltaMass>' e.g. \"-fm 'C:57.021464'\" would indicate a"
        " fixed carbamidomethylation of C."
        " Note: for fixed modifications on the same aminoacid, we only consider the first one."
        " For the N- and C-Terminus (for protein or peptide) use \"npepterm\", \"nproterm\" or \"cpepterm\", "
        "\"cproterm\""
    )
    group.add_argument(
        "--variable_mod", "-vm", type=_parse_mod, action="append",
        help="Apply a variable modification on a special aminoacid."
        " You can apply multiple variable modifications BUT only one modification per aminoacid is currently allowed."
        " The form should be '<AminoAcid>:<DeltaMass>' e.g. \"-vm 'M:15.994915'\" would indicate a variable oxidation"
        " of M."
        " For the N- and C-Terminus (for protein or peptide) use \"npepterm\", \"nproterm\" or \"cpepterm\", "
        "\"cproterm\""
    )

    # Flag to check if generated graphs are correctly generated
    group.add_argument(
        "--verify_graph", "--verify", default=False, action="store_true",
        help="Set this flag to perform a check whether the graph was generated correctly. Here we explicitly check "
        "for parallel edges, for DAG and other properties."
    )

    # Arguments for graph processing/digestion
    group.add_argument(
        "--digestion", "-d", type=str.lower, action="append",
        choices=["gluc", "trypsin", "skip", "full"],
        help="Set the digestion method. The full digestion cleaves at every edge, which can be useful for retrieving "
        "all possible peptides with arbitrary cutting points. The digestion method skip skips the digestion "
        "completely. You can use multiple digestions at once! Default: Trypsin"
    )
    group.add_argument(
        "--no_merge", "-nm", default=False, action="store_true",
        help="Set this flag to skip the merging process for chains of nodes and edges into a single node. Setting "
        "this option could drastically increase the size of the graph, especially its depth."
    )
    group.add_argument(
        "--no_collapsing_edges", default=False, action="store_true",
        help="Set this flag to indicate, that parallel edges should not be collapsed. This might be useful "
        "in some cases and might give insights if some errors occur. It can be especially used "
        "with the graph verification flag."
    )

    # Arguments for node and edge weights
    group.add_argument(
        "--annotate_mono_weights", "-amw", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the monoisotopic weights. (Values are taken from "
        "the mass dictionary)"
    )
    group.add_argument(
        "--annotate_avrg_weights", "-aaw", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the average weights. (Values are taken from "
        "the mass dictionary)"
    )
    group.add_argument(
        "--annotate_mono_weight_to_end", "-amwe", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the monoisotopic end weights. This weight informs about "
        "how much weight is at least left to get to the end Node. NOTE: Applying this, also sets the monoisotopic "
        "weights"
    )
    group.add_argument(
        "--annotate_avrg_weight_to_end", "-aawe", default=False, action="store_true",
        help="Set this to annotate nodes and edges with the average end weights. This weight informs about "
        "how much weight is at least left to get to the end Node. NOTE: Applying this, also sets the average weights"
    )
    group.add_argument(
        "--mass_dict_type", "-mdt",
        type=lambda s: int if s.lower() == "int" else (float if s.lower() == "float" else None), default="int",
        choices=[int, float], metavar="{int,float}",
        help="Set the type of the mass dictionary for amino acid. Default is set to int"
    )
    group.add_argument(
        "--mass_dict_factor", "-mdf", type=float, default=1000000000,
        help="Set the factor for the masses inside the mass_dictionary. The default is set to 1 000 000 000, "
        "so that each mass can be converted into integers."
    )
    group.add_argument(
        "--queue_size", "-qs", type=int, default=30000,
        help="Set the size of the queues ('reading of entries'-, 'writing of entries'- and statistics-queue), "
        "default is 30000"
    )


def add_statistics(group):
    group.add_argument(
        "--calc_num_possibilities", "-cnp", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This uses a dynamic programming approach to calculate this in an efficient manner."
    )
    group.add_argument(
        "--calc_num_possibilities_miscleavages", "-cnpm", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This returns a list, sorted by the number of miscleavages (beginning at 0). "
        "Example: Returns: [1, 3, 5, 2] -> 1 path with 0 miscleavages, 3 paths with 1 miscleavage, 5 paths "
        "with 2 miscleavages, etc ... This uses a dynamic programming approach to calculate this in an efficient "
        "manner. NOTE: This may get memory heavy, depending on the proteins (especially on Titin)"
    )
    group.add_argument(
        "--calc_num_possibilities_hops", "-cnph", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This returns a list, sorted by the number of hops (beginning at 0). "
        "Example: Returns: [0, 3, 5, 2] -> 0 paths with 0 hops, 3 paths with 1 hop, 5 paths "
        "with 2 hops, etc ... This uses a dynamic programming approach to calculate this in an efficient "
        "manner. NOTE: This mis even more memory heavy then binning on miscleavages. Of course it depends "
        "on the proteins (especially on Titin) NOTE: The dedicated start and end node is not counted here. "
        "If you traverse a graph, expect +2 more nodes in a path!"
    )

    group.add_argument(
        "--calc_num_possibilities_variant", "-cnpvar", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This returns a list, sorted by the number of variants (beginning at 0). "
        "Similar to misclavages"
    )

    group.add_argument(
        "--calc_num_possibilities_mutagen", "-cnpmut", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This returns a list, sorted by the number of mutagens (beginning at 0). "
        "Similar to misclavages"
    )

    group.add_argument(
        "--calc_num_possibilities_conflict", "-cnpcon", default=False, action="store_true",
        help="If this is set, the number of all possible (non repeating) paths from the start to the end node will"
        " be calculated. This returns a list, sorted by the number of conflicts (beginning at 0). "
        "Similar to misclavages"
    )

    # Add replace funcitonality to CLI
    def _list_to_func_map(input: str):
        """
        Choosing between possible or_count strategies for cnp
        """
        if "min" == input.lower():
            return min
        elif "max" == input.lower():
            return max
        else:
            return None
    group.add_argument(
        "--calc_num_possibilites_or_count", "-cnp_or_count", choices=[min, max],
        type=_list_to_func_map, action="store", default=min,
        help="Substitute amino acids in graphs by other amino acids. This could be useful to replace"
        " e.g. 'J' with 'I' and 'L'. This parameter can be then provided as: 'J->I,L'. Multiple replacements"
        " are allowed and are executed one after another. NOTE: only from ONE amino acid multiple amino acids can"
        " be substituted. So only the format: 'A->B[,C]*' is allowed!"
    )

    # TODO one parameter is missing for cnp con/mut/var


def add_graph_exports(group):
    group.add_argument(
        "--export_output_folder", "-eo", default=os.path.join(os.getcwd(), "exported_graphs"), type=str,
        help="Set the output folder to specify the dirctory of exported graphs (dot, graphml, gml) NOTE: It will "
        "overwrite exisiting files. Default is set the current working directory"
    )
    group.add_argument(
        "--export_in_directories", "-edirs", default=False, action="store_true",
        help="Set this flag to export files in directories (coded by accession) instead of directly by "
        "the accession name. This could be useful if millions of proteins are added into this tool, since "
        "a folder with millions of entries can be problematic in some cases."
    )
    group.add_argument(
        "--export_dot", "-edot", default=False, action="store_true",
        help="Set this flag to export a dot file for each protein"
    )
    group.add_argument(
        "--export_csv", "-ecsv", default=False, action="store_true",
        help="Set this flag to export a nodes-/edges-csv file for each protein"
    )
    group.add_argument(
        "--export_large_csv", "-elcsv", default=False, action="store_true",
        help="Set this flag to export a large nodes/edges-csv file, which contains every protein."
    )
    group.add_argument(
        "--export_graphml", "-egraphml", default=False, action="store_true",
        help="Set this flag to export a GraphML file for each protein. This is the recommended export method."
    )
    group.add_argument(
        "--export_gml", "-egml", default=False, action="store_true",
        help="Set this flag to export a GML file for each protein"
    )
    group.add_argument(
        "--export_pickle", "-epickle", default=False, action="store_true",
        help="Set this flag to export a Pickle file for each protein"
    )
    group.add_argument(
        "--export_pcsr", "-epcsr", default=False, action="store_true",
        help="Set this flag to export protein-graph with specific attributes in CSR-format"
    )
    group.add_argument(
        "--export_pcsr_pdb_entries", "-epcsr_pdbs", type=int, default=10,
        help="Set the number of pdb entries per node. Defaults to 10"
    )
    group.add_argument(
        "--export_large_pcsr", "-elpcsr", default=False, action="store_true",
        help="Set this flag to export protein-graph with specific attributes in CSR-format"
    )
    group.add_argument(
        "--export_large_pcsr_pdb_entries", "-elpcsr_pdbs", type=int, default=10,
        help="Set the number of pdb entries per node. Defaults to 10"
    )
    group.add_argument(
        "--export_binary_pcsr", "-ebpcsr", default=False, action="store_true",
        help="Set this flag to export protein-graph with specific attributes in a binary CSR-format"
    )
    group.add_argument(
        "--export_binary_pcsr_pdb_entries", "-ebpcsr_pdbs", type=int, default=10,
        help="Set the number of pdb entries per node. Defaults to 10"
    )
    group.add_argument(
        "--export_large_binary_pcsr", "-elbpcsr", default=False, action="store_true",
        help="Set this flag to export protein-graph with specific attributes in a binary CSR-format"
    )
    group.add_argument(
        "--export_large_binary_pcsr_pdb_entries", "-elbpcsr_pdbs", type=int, default=10,
        help="Set the number of pdb entries per node. Defaults to 10"
    )
    count_features = [
        "INIT_MET", "VARIANT", "VAR_SEQ", "SIGNAL", "MUTAGEN", "CONFLICT",
        "PEPTIDE", "PROPEP", "CHAIN", "VARMOD", "FIXMOD"
    ]
    group.add_argument(
        "--pcsr_feature_to_count", "-pcsr_ftc", choices=count_features, type=str.upper, action="append",
        help="Set the features which should be counted into the VC on a graph for all exportes in the pcsr. The values"
        " here correspond to the actual FT name provided in SwissProt-EMBL. Default is set to use only 'VARIANT',"
        " since most proteins are complicated due to many variants. Currently parsable features are: " +
        ", ".join(count_features) + ". Use it as follows to only select specific ones: '-pcsr_ftc VARIANT -pcsr_ftc"
        " MUTAGEN' to count mutagens and variants or '-pcsr_ftc PEPTIDE' to count only peptides in VC. This parameter"
        " can be used for later traversal to limit further exploration on specific nodes. E.G.: Only Paths from s to"
        " e, with maximum x many features (which can be checked on VC in pcsr)"
    )


def add_redis_graph_export(group):
    group.add_argument(
        "--export_redisgraph", "-eredisg", default=False, action="store_true",
        help="Set this flag to export to a redis-server having the module RedisGraph loaded."
    )
    group.add_argument(
        "--redisgraph_host", type=str, default="localhost",
        help="Set the host name for the redis-server having the module RedisGraph. Default: localhost"
    )
    group.add_argument(
        "--redisgraph_port", type=int, default=6379,
        help="Set the port for the redis-server having the module RedisGraph. Default: 6379"
    )
    group.add_argument(
        "--redisgraph_graph", type=str, default="proteins",
        help="Set the graph name on the redis-server having the module RedisGraph. Default 'proteins'"
    )


def add_postgres_graph_export(group):
    group.add_argument(
        "--export_postgres", "-epg", default=False, action="store_true",
        help="Set this flag to export to a postgresql server."
        "NOTE: This will try to create the tables 'nodes' and 'edges' on a specified database."
        " Make sure the database in which the data should be saved also exists."
    )
    group.add_argument(
        "--postgres_host", type=str, default="127.0.0.1",
        help="Set the host name for the postgresql server. Default: 127.0.0.1"
    )
    group.add_argument(
        "--postgres_port", type=int, default=5433,
        help="Set the port for the postgresql server. Default: 5433"
    )
    group.add_argument(
        "--postgres_user", type=str, default="postgres",
        help="Set the username for the postgresql server. Default: postgres"
    )
    group.add_argument(
        "--postgres_password", type=str, default="developer",
        help="Set the password for the postgresql server. Default: developer"
    )
    group.add_argument(
        "--postgres_database", type=str, default="proteins",
        help="Set the database which will be used for the postgresql server. Default: proteins"
    )


def add_mysql_graph_export(group):
    group.add_argument(
        "--export_mysql", "-emysql", default=False, action="store_true",
        help="Set this flag to export to a MySQL server."
        "NOTE: This will try to create the tables 'nodes' and 'edges' on a specified database."
        " Make sure the database in which the data should be saved also exists."
    )
    group.add_argument(
        "--mysql_host", type=str, default="127.0.0.1",
        help="Set the host name for the MySQL server. Default: 127.0.0.1"
    )
    group.add_argument(
        "--mysql_port", type=int, default=3306,
        help="Set the port for the MySQL server. Default: 3306"
    )
    group.add_argument(
        "--mysql_user", type=str, default="root",
        help="Set the username for the MySQL server. Default: root"
    )
    group.add_argument(
        "--mysql_password", type=str, default="",
        help="Set the password for the MySQL server. Default: <empty>"
    )
    group.add_argument(
        "--mysql_database", type=str, default="proteins",
        help="Set the database which will be used for the MySQL server. Default: proteins"
    )


def add_cassandra_export(group):
    group.add_argument(
        "--export_cassandra", "-ecassandra", default=False, action="store_true",
        help="Set this flag to export to a Cassandra server."
    )
    group.add_argument(
        "--cassandra_host", default="127.0.0.1", type=str,
        help="Set the address/host of a cassandra server (One server from the cluster is sufficient). "
        "Default is set to 127.0.0.1"
    )
    group.add_argument(
        "--cassandra_port", default=9042, type=int,
        help="Set the port which should be used to connect to the cassandra host. "
        "Default is set to 9042"
    )
    group.add_argument(
        "--cassandra_keyspace", default="graph", type=str,
        help="Set the keyspace where ProtGraph can operate on. If the keyspace does not exist. ProtGraph will attempt "
        " to create it.  Default keyspace is 'graph'"
    )
    group.add_argument(
        "--cassandra_chunk_size", default=50, type=int,
        help="Set the size of batches, which should be then sent to cassandra. "
        " The default is set specifically low, since, cassandra has configured a low batch size (50kb)"
    )


def add_peptide_export_traversal_options(group):
    group.add_argument(
        "--pep_hops", type=int, default=-1,
        help="Set the number of hops (max length of path) which should be used to get paths "
        "from a graph. NOTE: the larger the number the more memory may be needed. This depends on "
        "the protein which currently is processed. Default is set to '-1', so all lengths are considered."
    )
    group.add_argument(
        "--pep_miscleavages", type=int, default=-1,
        help="Set this number to filter the generated paths by their miscleavages."
        "The protein graphs do contain infomration about 'infinite' miscleavages and therefor also return "
        "those paths/peptides. If setting (default) to '-1', all results are considered. However you can limit the "
        "number of miscleavages, if needed."
    )
    group.add_argument(
        "--pep_skip_x",  default=False, action="store_true",
        help="Set this flag to skip to skip all entries, which contain an X"
    )
    group.add_argument(
        "--pep_use_igraph",  default=False, action="store_true",
        help="Set this flag to use igraph instead of netx. "
        "NOTE: If setting this flag, the peptide generation will be considerably faster "
        "but also consumes much more memory. Also, the igraph implementation DOES NOT go "
        "over each single edge, so some (repeating results) may never be discovered when using "
        "this flag."  # TODO see issue: https://github.com/igraph/python-igraph/issues/366
    )
    group.add_argument(
        "--pep_min_pep_length", type=int, default=0,
        help="Set the minimum peptide length to filter out smaller existing path/peptides. "
        "Here, the actual number of aminoacid for a peptide is referenced. Default: 0"
    )
    group.add_argument(
        "--pep_batch_size", type=int, default=25000,
        help="Set the batch size. This defines how many peptides are inserted at once. "
        "Default: 25000"
    )
    group.add_argument(
        "--pep_min_weight", type=float, default=-1,
        help="Set the minimum peptide weight in Da which should be considered. Set negative to ignore. "
        "Default: -1"
    )
    group.add_argument(
        "--pep_max_weight", type=float, default=-1,
        help="Set the maximum peptide weight in Da which should be considered. Set negative to ignore. "
        "Default: -1"
    )


def add_postgres_peptide_export(group):
    group.add_argument(
        "--export_peptide_postgres", "-epeppg", default=False, action="store_true",
        help="Set this flag to export peptides (specifically paths) to a postgresql server."
        "NOTE: This will try to create the tables 'accessions' and 'peptides' on a specified database."
        " Make sure the database in which the data should be saved also exists. If problems occur, try "
        "to delete the generated tables and retry again."
    )
    group.add_argument(
        "--pep_postgres_host", type=str, default="127.0.0.1",
        help="Set the host name for the postgresql server. Default: 127.0.0.1"
    )
    group.add_argument(
        "--pep_postgres_port", type=int, default=5433,
        help="Set the port for the postgresql server. Default: 5433"
    )
    group.add_argument(
        "--pep_postgres_user", type=str, default="postgres",
        help="Set the username for the postgresql server. Default: postgres"
    )
    group.add_argument(
        "--pep_postgres_password", type=str, default="developer",
        help="Set the password for the postgresql server. Default: developer"
    )
    group.add_argument(
        "--pep_postgres_database", type=str, default="proteins",
        help="Set the database which will be used for the postgresql server. Default: proteins"
    )
    group.add_argument(
        "--pep_postgres_no_duplicates",  default=False, action="store_true",
        help="Set this flag to not insert duplicates into the database. "
        "NOTE: Setting this value decreases the performance drastically"
    )


def add_sqlite_peptide_export(group):
    # TODO
    group.add_argument(
        "--export_peptide_sqlite", "-epepsqlite", default=False, action="store_true",
        help="Set this flag to export peptides (specifically paths) to a sqlite server."
        "NOTE: This will try to create the tables 'accessions' and 'peptides' on a specified database."
        " Make sure the database in which the data should be saved also exists. If problems occur, try "
        "to delete the generated tables and retry again."
    )
    group.add_argument(
        "--pep_sqlite_database", type=str, default="peptides.db",
        help="Set the database file which will be used for the sqlite server. Default: proteins"
    )
    group.add_argument(
        "--pep_sqlite_output_dir", type=str, default=os.getcwd(),
        help="Set the output directory which will be used for the sqlite file. Defaults to ./"
    )
    group.add_argument(
        "--pep_sqlite_non_compact",  default=False, action="store_true",
        help="Set this flag to to have a non compact table repr"
    )
    group.add_argument(
        "--pep_sqlite_compression_level", type=int, default=6,
        help="Set the compression level of the blob (zlib). Note: This value "
        "is only used of 'use_blobs' is set. Default: 6"
    )
    group.add_argument(
        "--pep_sqlite_use_blobs",  default=False, action="store_true",
        help="Set this flag to save binary compressed blobs instead of strings in the database (reducing its size)."
    )


def add_mysql_peptide_export(group):
    group.add_argument(
        "--export_peptide_mysql", "-epepmysql", default=False, action="store_true",
        help="Set this flag to export peptides (specifically paths) to a MySQL server."
        "NOTE: This will try to create the tables 'accessions' and 'peptides' on a specified database."
        " Make sure the database in which the data should be saved also exists. If problems occur, try "
        "to delete the generated tables and retry again."
    )
    group.add_argument(
        "--pep_mysql_host", type=str, default="127.0.0.1",
        help="Set the host name for the mysql server. Default: 127.0.0.1"
    )
    group.add_argument(
        "--pep_mysql_port", type=int, default=3306,
        help="Set the port for the mysql server. Default: 3306"
    )
    group.add_argument(
        "--pep_mysql_user", type=str, default="root",
        help="Set the username for the mysql server. Default: root"
    )
    group.add_argument(
        "--pep_mysql_password", type=str, default="",
        help="Set the password for the mysql server. Default: ''"
    )
    group.add_argument(
        "--pep_mysql_database", type=str, default="proteins",
        help="Set the database which will be used for the mysql server. Default: proteins"
    )
    group.add_argument(
        "--pep_mysql_no_duplicates",  default=False, action="store_true",
        help="Set this flag to not insert duplicates into the database. "
        "NOTE: Setting this value decreases the performance drastically"
    )


def add_fasta_peptide_export(group):
    group.add_argument(
        "--export_peptide_fasta", "-epepfasta", default=False, action="store_true",
        help="Set this flag to export peptides into a single fasta file."
    )
    group.add_argument(
        "--pep_fasta_out", default=os.path.join(os.getcwd(), "peptides.fasta"),
        type=str,
        help="Set the output file for the peptide fasta export. Default: '${pwd}/peptides.fasta'. "
        "NOTE: This will overwrite existing files."
    )


def add_trie_peptide_export(group):
    group.add_argument(
        "--export_peptide_trie", "-epeptrie", default=False, action="store_true",
        help="Set this flag to export peptides into a single fasta file. "
        "NOTE: This exports peptides in a trie sturcure on the filesystem. Make sure your "
        "to be exported filesystem supports generating millions of folders and files! (preferably XFS)"
    )
    group.add_argument(
        "--pep_trie_folder_out", default=os.path.join(os.getcwd(), "exported_peptides"),
        type=str,
        help="Set the output folder for the peptide export. Default: '${pwd}/peptides.fasta'. "
        "NOTE: This will append existing files."
    )


def add_citus_peptide_export(group):
    group.add_argument(
        "--export_peptide_citus", "-epepcit", default=False, action="store_true",
        help="Set this flag to export peptides (specifically paths) to a postgresql-citus server."
        "NOTE: This will try to create the tables 'accessions' and 'peptides' on a specified database."
        " Make sure the database in which the data should be saved also exists. If problems occur, try "
        "to delete the generated tables and retry again."
    )
    group.add_argument(
        "--pep_citus_host", type=str, default="127.0.0.1",
        help="Set the host name for the postgresql server with citus. Default: 127.0.0.1"
    )
    group.add_argument(
        "--pep_citus_port", type=int, default=5432,
        help="Set the port for the postgresql server with citus. Default: 5432"
    )
    group.add_argument(
        "--pep_citus_user", type=str, default="postgres",
        help="Set the username for the postgresql server with citus. Default: postgres"
    )
    group.add_argument(
        "--pep_citus_password", type=str, default="developer",
        help="Set the password for the postgresql server with citus. Default: developer"
    )
    group.add_argument(
        "--pep_citus_database", type=str, default="proteins",
        help="Set the database which will be used for the postgresql server with citus. Default: proteins"
    )
    group.add_argument(
        "--pep_citus_no_duplicates",  default=False, action="store_true",
        help="Set this flag to not insert duplicates into the database. "
        "NOTE: Setting this value decreases the performance drastically"
    )


def add_gremlin_graph_export(group):
    group.add_argument(
        "--export_gremlin", "-egremlin", default=False, action="store_true",
        help="Set this flag to export the graphs via gremlin to a gremlin server."
        "NOTE: The export is very slow, since it executes each node as a single query "
        "(tested on JanusGraph and Apache Gremlin Server). This exporter is not well implemented and may not work. "
        "This is due to difficulties implementing such an exporter in a global manner. "
        "To reduce the number of errors: Try to have a stable connection to the gremlin-server and also allocate "
        "enough resources for it, so that it can process the queries quick enough."
    )
    group.add_argument(
        "--gremlin_url", type=str, default="ws://localhost:8182/gremlin",
        help="Set the url to the gremlin URL (no authentication). Default: 'ws://localhost:8182/gremlin'"
    )
    group.add_argument(
        "--gremlin_traversal_source", type=str, default="g",
        help="Set the traversal source for remote. Default 'g'"
    )
