# [WIP] ProtGraph - A Graph-Generator for Proteins

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
FT                   /note="L -> K (in Example)"
FT                   /evidence="EXAMPLE Variant 3"
FT                   /id="VAR_XXXXX3"
SQ   SEQUENCE   8 AA;  988 MW;  XXXXXXXXXXXXXXXX CRC64;
     MPROTEIN                   
```

This entry contains the canonical sequence (described in section `SQ`), which is usually exported in fasta files and then used by other tools. But such entries provide more information, which are usually not considered. Especially in this example, a `SIGNAL` peptide is present, which may lead to a peptide, which is not covered when cleaving (e.g. via the digestion enzyme Trypsin).
Additionally, 3 `VARIANT`s are present in this entry, informing about substituted aminoacids, which again may yield peptides that are not covered by the canonical sequence or by digestion.

With `ProtGraph`, a graph can be generated. In this example, it would look like follows:
![Example Graph](resources/example_graph.png)

With such a graph representation, the earlier mentioned peptides, which might not be covered are now covered and could be retrieved, by traversing from the dedicated `s`tart node to the dedicated `e`nd node. (These start and end nodes are internally denoted by: `__start__` and `__end__`).
The `black drawn` edges (with the label `TRYPSIN`) indicate that this graph was additionally digested via the enzyme Trypsin (default) leading to a compact directed and acyclic graph (in short DAG). From such DAGs, all possible peptides and further statistics can be retrieved.

ProtGraph will generate a statistics file when generating a graph from an UniProt-entry, containing various information about the protein/entry which were retrieved during generation (partially on the fly).

Here are some examples, which can be retrieved from the statistics:

If executing via the flag `-cnp`, the total number of possible non-repeating paths between `s` and `e` are calculated. In this example, `46` paths are possible. With the number of possible paths, we know that the graph contains `46` peptides (some of them repeating).
It may be interesting to know how many (repeating) peptides with a certain number of miscleavages are present in a protein. To calculate this statistic, simply provide the flag `-cnpm`. In this example: `23` peptides with 0 miscleavages, `19` with 1 miscleavage and `4` with 2 miscleavages are present.

By combining multiple flags, it is possible to retrieve the number of  peptides having a specific length (counting aminoacids). With the flags `-cnph` and `-nm`, ProtGraph would write into the statistics file the needed information to retrieve the following:
| Peptide Length | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|----------------|---|---|---|---|---|---|---|---|
| #Peptides      | 3 | 5 | 8 | 8 | 6 | 4 | 8 | 4 |


Dividing this table by `46` would give us the distribution of the peptide lengths for the given protein.



## Setting up ProtGraph

First clone this project and enter its directory via:

```shell
git clone git@github.com:mpc-bioinformatics/ProtGraph.git
cd ProtGraph/
```

Depending on how you want to use this project follow either the usage instructions for only using Protgraph or the development and tests instructions if you want to use it in a python project or even develop something in ProtGraph itself.

### Installing via pip (for usage)

If you only want to use ProtGraph from the command line, it is sufficient to simply execute: `pip install --user .`

The "binary" `protgraph` should be available in your command line after it finished. You can now convert UniProt-entries into graphs!

#### Troubleshooting

ProtGraph has many dependencies which are mostly due to the export functionalities. For the dependency `psycopg2`, it may be neccessary to install PostgreSQL on your operating system, so that the wheel can be built. It should be sufficient to install it as follows `sudo pacman -S postgres` (adapt this for `apt` and alike).

If the command `protgraph` cannot be found, make sure that you have included the python executables folder in your environment variable `PATH`. (If pip was executed with the flag `--user`, then the binaries should be available in the folder: `~/.local/bin`)

### Installing via pipenv and pyenv (for development and tests)

To set up ProtGraph make sure you have `pipenv` and `pyenv` installed. (via `pacman` `apt` and alike)


The dependencies (and possibly Python) need to be installed first. This can be done with:
> pipenv install

Afterwards execute to install the module ProtGraph:
> pipenv run pip install -e .

If both commands finished succesfully, activate the environment via:
> pipenv shell

You should be able to run `protgraph` now. Via the flag `-e` the code from the projects directory is used, so you can edit the files directly and test the changes via your debugger of choice or via the command `protgraph`.


#### Testing

Make sure you are at the root folder of this project. Then execute `pytest .` to run the tests for ProtGraph.

## Usage

To see an overview of possible parameters and options use: `protgraph --help`. The `help`-output contains brief descriptions for each flag.

### Example calls

Let's use the provided example `e_coli.dat` located in `examples/e_coli.dat` (Or any other `.txt` or `.dat` file).

The graph generation can be executed via: `protgraph examples/e_coli.dat`. This will generate the graphs and the additional statistics file. You can inspect the statistics file after it has finished. This should only take a few seconds/minutes.

---

During execution it could be seen that now information was provided, when ProtGraph finishes.

Sadly, `Biopython` does not provide information of how many entries are available in a `.txt` or `.dat` file. Therefore this information needs to be provided via another parameter: `protgraph --num_of_entries 9434 examples/e_coli.dat` (you can also use `-n`).

To retrieve the number of entries beforehand, you could e.g. use `cat examples/e_coli.dat | grep "^//" | wc -l`. It is also possible to add multiple files. The number of entries for each file then need to be summed: `protgraph -n 18868 examples/e_coli.dat examples/e_coli.dat`

---

If to many (or to few) processes are executed, then it can be adjusted via the parameter `--num_of_processes` or `-np`. E.G. `protgraph --num_of_processes 3 --num_of_entries 9434 examples/e_coli.dat` will use 4 (`3` + 1 reading process) processes.

---

To fully annotate the graphs with weights and to retrieve currently all available statistics, use the following: `protgraph -amwe -aawe -cnp -cnpm -cnph -n 9434 examples/e_coli.dat`

## Exporting graphs

While executing ProtGraph, generated graphs are not saved.

This is the default behaviour of `protgraph`. It excludes the generated graphs, since those can explode in size and the disk space on machines  may also be limited. Currently a few export functionalities are available and it is planned to extend this functionality.

### File Exports

To export each protein graph into a file, simply set the flags `-edot`, `-egraphml` and/or `-egml`, which will create the corresponding dot, GraphML or GML files into a output folder.
Exporting to GraphML is recommended since this is the only export method able to serialize all properties in a file which can be set in the graph.

### Database Exporters

The database exporters are currently under development. Two exports are currently available. `protgraph` allows to export the generated graphs into PostgreSQL as well as into RedisGraph. For configuring such export, please look into the `--help` output.

Furthermore an export (experimental) via gremlin is provided as well as an peptide export into PostgreSQL (NOTE: Peptide exports may generate unmanagable amounts of peptides).

Note: In PostgreSQL a database should be created with no tables to ensure that ProtGraph can generate the tables itself without errors.

Note: Currently no peptide export is available, but it is planned to be added into ProtGraph. However, using the peptide export to PostgreSQL, some can easily script a peptide export from the database and the corresponding graph files by iterating over each entry.


## Missing Functionalities

* Other export functionalities
* Peptide export scripts/possibility for ProtGraph
