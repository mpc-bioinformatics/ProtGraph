# ProtGraph - A Protein-Graph-Generator

[![install with pip](https://img.shields.io/badge/Install%20with-pip-brightgreen)](https://pypi.org/project/protgraph/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/protgraph?label=PyPI%20-%20Downloads)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/protgraph/README.html)
![Conda](https://img.shields.io/conda/dn/bioconda/protgraph?label=BioConda%20-%20Downloads)

## Summary

ProtGraph in short is a python-package, which allows to convert protein-entries from the UniProtKB to so-called protein-graphs. We use the SP-EMBL-Entries, provided by UniProtKB via `*.txt` or `*.dat`-files, and parse the available feature-information. In contrast to a FASTA-file-entry of a protein, a SP-EMBL-file-entry is more detailed. SP-EMBL-files not only contain the canonical sequence, but also contain isoform-sequences, specifically cleaved peptides, like SIGNAL-peptides, PROPEPtides (or simply: PEPTIDEs, e.g. ABeta 40/42), variationally altered aminoacid-sequences (VARIANTs, mutations or sequence-CONFLICTs) and even more. The SP-EMBL-files are parsed, using biopyton, and corresponding protein-graphs, containing the additional feature-information are generated with igraph. The resulting graphs contain all possibly resulting protein-sequences from the canonical as well as from the feature-information and can also be further processed (e.g. digestion, to contain peptides).

So, what can we do with ProtGraph? First, the protein-graphs can be saved by ProtGraph in various formats. These can be then opened via python/R/C++/... or any other tool (e.g. Gephi) and algorithms, statistics, visualizations, etc... can be applied. ProtGraph then acts solely as a converter of protein-entries to a graph-structure. Second, while generating protein-graphs, a statistics-file is generated, containing various on-the-fly retrievable information about the protein-graph. We can calculate the number of nodes/edges/features within a protein-graph, the number of protein-/peptide-sequences contained and even binned by specific attributes. These can give a quick overview e.g. of the tryptic search space while considering all feature-information of the provided E.Coli-Proteome (containing in total XX peptides). Lastly, the protein-graphs per se are not usefull, especially in identification. Therefore, we extended ProtGraph to optionally convert protein-graph into FASTA-entries, to be used for identification, to enable searches of feature-induced-peptides.

Curious what we do with ProtGraph and its output? Visit materials_and_posters for an explanation of the protein-graph-generation and further materials. Further below in the README.md, additional examples are provided.

## Generating a Protein-Graph from an Example-Entry and retrieving Statistics

This project aims to efficiently generate graphs of proteins using text entries from [UniProt](https://www.uniprot.org/). These text files usually end with the extension: `.dat` or `.txt`.

We utilize the tools: [BioPython](https://biopython.org/) and [igraph](https://igraph.org/python/) to generate so called directed and acyclic graphs (DAGs) of proteins.

In essence, ProtGraph is doing the following:

Lets assume that an entry has been exported as text from UniProt, containing the following content:

```txt
ID   MY_CUSTOM_PROTEIN             Reviewed;         8 AA.
AC   QXXXXX

...

FT   INIT_MET        1
FT                   /note="Removed"
FT                   /evidence="EXAMPLE Initiator Methionine"
FT   SIGNAL          1..4
FT                   /evidence="EXAMPLE Signal Peptide"
FT   VARIANT         3
FT                   /note="Missing (in Example)"
FT                   /evidence="EXAMPLE Variant 1"
FT                   /id="VAR_XXXXX1"
FT   VARIANT         4
FT                   /note="O -> L (in Example)"
FT                   /evidence="EXAMPLE Variant 2"
FT                   /id="VAR_XXXXX2"
FT   VARIANT         5
FT                   /note="T -> K (in Example)"
FT                   /evidence="EXAMPLE Variant 3"
FT                   /id="VAR_XXXXX3"
SQ   SEQUENCE   8 AA;  988 MW;  XXXXXXXXXXXXXXXX CRC64;
     MPROTEIN                   
```

This entry contains the canonical sequence (described in section `SQ`), which is usually exported in fasta files and then used by other tools. But such entries provide more information, which are usually not considered. Especially in this example, a `SIGNAL` peptide is present, which may lead to a peptide, which is not covered when cleaving (e.g. via the digestion enzyme Trypsin).
Additionally, 3 `VARIANT`s are present in this entry, informing about substituted aminoacids, which again may yield peptides that are not covered by the canonical sequence or by digestion.

With `ProtGraph`, a graph can be generated. In this example, it would look like follows:
![Example Graph](https://raw.githubusercontent.com/mpc-bioinformatics/ProtGraph/master/resources/example_graph.png)

With such a graph representation, the earlier mentioned peptides, which might not be covered are now covered and could be retrieved, by traversing from the dedicated `s`tart node to the dedicated `e`nd node. (These start and end nodes are internally denoted by: `__start__` and `__end__`).
The `orange drawn` edges (with the label `TRYPSIN`) indicate that this graph was additionally digested via the enzyme Trypsin (default) leading to a compact directed and acyclic graph (in short DAG). From such DAGs, all possible peptides and further statistics can be retrieved.

ProtGraph will generate a statistics file when generating a graph from an UniProt-entry, containing various information about the protein/entry which were retrieved during generation (partially on the fly).

Here are some examples, which can be retrieved from the statistics:

If executing via the flag `-cnp`, the total number of possible non-repeating paths between `s` and `e` are calculated. In this example, `46` paths are possible. With the number of possible paths, we know that the graph contains `46` peptides (some of them possibly repeating).
It may be interesting to know how many (repeating) peptides with a certain number of miscleavages are present in a protein. To calculate this statistic, simply provide the flag `-cnpm`. In this example: `23` peptides with 0 miscleavages, `19` with 1 miscleavage and `4` with 2 miscleavages are present.

By combining multiple flags, it is possible to retrieve the number of  peptides having a specific length (counting aminoacids). With the flags `-cnph` and `-nm`, ProtGraph would write into the statistics file the needed information to retrieve the following:

| Peptide Length | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|----------------|---|---|---|---|---|---|---|---|
| #Peptides      | 3 | 5 | 8 | 8 | 6 | 4 | 8 | 4 |

Dividing this table by `46` would give us the distribution of the peptide lengths for the given protein.

## Setting up ProtGraph

You can download Protgraph directly from [pypi](https://pypi.org/project/protgraph).

Simply execute: `pip install protgraph`  or `pip install --user protgraph`

It is also avialable in [bioconda](https://bioconda.github.io/). The installation instruction can be found here: [protgraph](https://anaconda.org/bioconda/protgraph)

### Installing via pip from source (for usage and development)

First clone this project and enter its directory via:

```shell
git clone git@github.com:mpc-bioinformatics/ProtGraph.git
cd ProtGraph/
```

Depending on how you want to use this project follow either the usage instructions or the development and tests instructions.

#### Installing for Usage

If you only want to use ProtGraph from the command line, it is sufficient to simply execute: `pip install --user .`. Alternatively you can execute `python setup.py install`.

The "binary" `protgraph` should be available in your command line after it finished. You can now convert UniProt-entries into graphs!

#### Installing for Development and Test

To set up ProtGraph make sure you have `pipenv` and `pyenv` installed. (via `pacman` `apt` and alike)

First install needed dependencies with:
> pipenv install && pipenv install -d

Afterwards execute the following to install the module ProtGraph inside the environment:
> pipenv run pip install -e .

If both commands finished succesfully, activate the environment via:
> pipenv shell

You should be able to run `protgraph` now. Via the flag `-e` the code from the projects directory is used, so you can edit the files directly and test the changes via your debugger of choice or via the command `protgraph`.

Tests can be executed as follows: Make sure you are at the root folder of this project. Then execute `pytest .` to run the tests for ProtGraph.

#### Troubleshooting

ProtGraph has many dependencies which are mostly due to the export functionalities. For the dependency `psycopg`, it may be neccessary to install PostgreSQL on your operating system, so that the wheel can be built. It should be sufficient to install it as follows `sudo pacman -S postgres` (adapt this for `apt` and alike).

If the command `protgraph` cannot be found, make sure that you have included the python executables folder in your environment variable `PATH`. (If pip was executed with the flag `--user`, then the binaries should be available in the folder: `~/.local/bin`)

## Usage

To see an overview of possible parameters and options use: `protgraph --help`. The `help`-output contains brief descriptions for each flag.

### Example calls

Let's use the provided example `e_coli.dat` located in `examples/e_coli.dat` (Or any other `.txt` or `.dat` file).

The graph generation can be executed via: `protgraph examples/e_coli.dat`. This will generate the graphs and the additional statistics file. You can inspect the statistics file after it has finished. This should only take a few seconds/minutes.

A progressbar was not shown during execution which is due to `Biopython` not providing information of how many entries are available in a `.txt` or `.dat` file. Therefore this information needs to be provided via another parameter: `protgraph --num_of_entries 9434 examples/e_coli.dat` (you can also use `-n`).

To retrieve the number of entries beforehand, you could e.g. use `cat examples/e_coli.dat | grep "^//" | wc -l`. It is also possible to add multiple files into ProtGraph. The number of entries for each file then needs to be summed: `protgraph -n 18868 examples/e_coli.dat examples/e_coli.dat`

If too many (or too few) processes are executed, then it can be adjusted via the parameter `--num_of_processes` or `-np`. E.G. `protgraph --num_of_processes 3 --num_of_entries 9434 examples/e_coli.dat` will use 4 (`3` + `1` reading) processes.

To fully annotate the graphs with weights and to retrieve currently all available statistics, use the following: `protgraph -amwe -aawe -cnp -cnpm -cnph -n 9434 examples/e_coli.dat`

---

Substitutions examples:

Protgraph is able to substitute amino acids if needed. E.G. it could be interesting to have the amino acids `I` and `L` instead of `J`.
Substitution can be done with the argument `--replace_aa` (in short `-raa`):
> protgraph -n -raa "J->I,L" examples/e_coli.dat

Protgraph is not limited to only substitute `J`. These substitutions are configurable by the user as long as these conform to the following format: `X->A[,B]*`. Multiple substitutions can be provided:
> protgraph -n -raa "J->I,L" -raa "I -> L" examples/e_coli.dat

Note: Substitutions are executed one after another which leads in this example to effectively substitues `J->L` and `I->L`.

In the following all common substitutions are listed
> protgraph -n -raa "B->D,N" examples/e_coli.dat
>
> protgraph -n -raa "J->I,L" examples/e_coli.dat
>
> protgraph -n -raa "Z->Q,E" examples/e_coli.dat
>
> protgraph -n -raa "X->A,C,D,E,F,G,H,I,K,L,M,N,O,P,Q,R,S,T,U,V,W,Y" examples/e_coli.dat

---

Fasta export examples:

Protgraph can be used to get all possible peptides from isoforms and canonoical for `e_coli.dat` into a fasta file:
> protgraph -n 9434 -ft NONE -epepfasta examples/e_coli.dat

Instead of only containing peptides, ProtGraph can get all possible proteins with isoforms and canonoical for `e_coli.dat` into a fasta file:
> protgraph -n 9434 -ft VAR_SEQ -d skip -epepfasta examples/e_coli.dat

It could also be interesing to get a fasta file of proteins with all variants, isoforms, cleaved or not cleaved signal peptide/initiator methionine:
> protgraph -n 9434 -ft VARIANT -ft VAR_SEQ -ft SIGNAL -ft INIT_MET -d skip -epepfasta examples/e_coli.dat

(This generates a 5.2G Fasta file!)

It could also be interesing to get a fasta file of proteins with all available features ProtGraph can parse:
> protgraph -n 9434 -ft ALL -d skip -epepfasta examples/e_coli.dat

or

> protgraph -n 9434 -d skip -epepfasta examples/e_coli.dat

Maybe it is neccessary to retrieve all possible peptide combinations of arbitrary cut peptides in proteins. This can be achieved via `-d full`. However, this would take a lot of time and a lot of disk space for the whole dataset (even for `e_coli.dat`). Instead we extract all possible peptides which do not have more than `100` aminoacids:
> protgraph -n 9434 -ft NONE -d full -epepfasta --pep_fasta_hops 100 examples/e_coli.dat

(This generates a 20G Fasta file with 244 985 079 peptides!)

Here the example with all possible peptide combinations for `e_coli.dat` to illustrate how quickly it can explode:
> protgraph -n 9434 -ft NONE -d full -epepfasta - examples/e_coli.dat

(This generates a 219G Fasta file with 694 844 270 peptides!)

It is advisable to do a dry run with the flags `-cnp`, `-cnpm` or `cnph` (or all of them) before exporting peptides/proteins to have an overview of how many peptides/proteins are currently represented by all graphs and if it is feasible to generate a Fasta file.

## ProtGraph Exporters

While executing ProtGraph, generated graphs are not saved. This is also true for peptides or proteins, which are represented by the graphs.

This is the default behaviour of `protgraph`. It excludes the generated graphs/proteins/peptides, since those can explode in size and the disk space on machines may also be limited. This is illustrated in the above examples.

However, those exporting functionalities can be used. Those differ in two categories: `-e*` for graph exports and `-epep*` for peptide/protein exports. For the peptide/protein exports limitations can be set on the graph-traversal itself.

**NOTE:** Graphs can contain unmanagable amounts of peptides/proteins. Do a dry run WITHOUT the export functionality first and examine the statistics output. Each peptide/protein exporter has multiple parameters to further limit the number of results.
Without a dry run it may happen that a protein like P04637 (P53 Human) with all possible peptides and variants is exported, which will very likely take up all disk space.

It is planned to add further export functionalities if needed.

### File Exports

To export each protein graph into a file, simply set the flags `-edot`, `-egraphml` and/or `-egml`, which will create the corresponding dot, GraphML or GML files into a output folder. With the flag `-epickle` it is also possible to generate a binary pickle file of the graph (which can be used by other python programs with `igraph`).
If processing many proteins (over 1 000 000) it is recommended to set the flag `-edirs` to save the graphs into multiple folders, since the underlying filesystem may be overwhelmed with too many files in one folder.

Exporting to GraphML is recommended since this is the only export method able to serialize all properties in a file which may be set in a graph. However, keep in mind that this export needs the most disk space out of the other graph exporters.

There is also the possibility to directly export graphs into a `fasta` file. This can be achieved via: `-epepfasta`.

### Database Exporters

ProtGraph offers multiple database-/storage-exporters. It is possible to export graphs into PostgreSQL/MySQL/RedisGraph and Gremlin (experimental). For PostgreSQL as well as MySQL two tables (`nodes`, `edges`) are generated and corresponding entries for each graph are added.

The export to RedisGraph, generates a graph which can be then accessed by other tools. The experimental export to Gremlin is in a working state, but is not further developed. It works with JanusGraph as well as with TinkerPop (Gremlin Server).

Additionally peptides/proteins can be exported via PostgreSQL as well as with MySQL. Both contain the option to skip duplicated entries, which can reduce the size of those databases.

For more information of exporters and other functionalities of ProtGraph, execute: `protgraph --help`.
