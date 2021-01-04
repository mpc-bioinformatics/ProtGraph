# WIP Graph-Generator for Proteins

This project aims to efficiently generate graphs of proteins from [UniProt](https://www.uniprot.org/) `.dat` or `.txt` files.

It does so by utilizing Biopython and igraph to generate directed acyclic graphs of Proteins 


[TODO IMAGE describing what happens from entry to graph]



## Setting up

Make sure you have `pipenv` and `pyenv` installed. (via `pacman` `apt` and alike)

First clone this project and enter its directory via: 
```shell
git clone <project-url>
cd <project-url>
```

Now the dependencies (and possibly Python) need to be installed. This can be done with:
> pipenv install

After everything is set up activate the environment via:
> pipenv shell

You are now ready to use this project!


## Usage
Currently, the entrypoint of this script is `main.py`. However, this may change in future.


To see an overview of possible parameters and options use: `python main.py --help`

---

Lets assume that `e_coli.txt` has been downloaded from Uniprot (as `.txt`)

The Graphgeneration can be executed via: `python main.py e_coli.txt`. However, we do not see how long it takes to finish. We only see the number of processed proteins per second. This tool can provide an estimation, when providing the number of entries. Currently (as for 04.01.2021), the E-Coli dataset has 9434 entries. Executing `python main.py --num_of_entries 9434 e_coli.txt` will print this (you can also use `-n`).


It is also possible to add multiple files. The number of entries for each file then need to be summed: `python main.py -n 18868 e_coli.txt e_coli.txt`


If to many (or to few) processes are executed, then it can be adjusted via `--num_of_processes` or `-np`. E.G. `python main.py --num_of_processes 3 --num_of_entries 9434 e_coli.txt ` will use 4 (3 + 1 reading process) processes.


To execute everything there can be, use the following: `python main.py -amew -aaew -cnp -n 9434 r_coli.dat`


Proteins are digested by Trypsin (default). This can be changed by setting the `--digestion` to something else. See help information for more.


--- 

A statistics file will be generated during each run. This statistics file contains various information, which were gathered on the fly or computed (if specified). 

## Missing Functionalities

* Export to graph files in folder
* Export to Postgres
* Other Export functionalities!





