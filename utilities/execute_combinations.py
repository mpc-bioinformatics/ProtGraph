import itertools

import protgraph

# This is a dummy file for iterating over each combination of specified features
# and getting the corresponding results in an automatic manner.

if __name__ == "__main__":

    used_features = [
        # (identifier, configuration_parameter)
        ("a", ["-ft", "INIT_MET"]),
        ("b", ["-ft", "SIGNAL"]),
        ("c", ["-ft", "VAR_SEQ"]),
        ("d", ["-ft", "VARIANT"]),
        ("e", ["-ft", "MUTAGEN"]),
        ("f", ["-ft", "CONFLICT"]),
        ("g", ["-ft", "NONE", "-raa", "B->D,N"]),
        ("h", ["-ft", "NONE", "-raa", "J->I,L"]),
        ("i", ["-ft", "NONE", "-raa", "Z->Q,E"])
    ]

    input_files = ["../examples/e_coli.dat"]
    statistics = ["-cnp", "-n", "9434", "-no_desc"]

    counter = 0
    for r in range(0, len(used_features) + 1):
        for c in itertools.combinations(used_features, r):
            counter += 1

            print("\nExecuting Iteration {} out of {}".format(counter, 2**(len(used_features))))
            if len(c) != 0:
                configuration = [y for x in c for y in x[1]]
            else:
                configuration = ["-ft", "NONE"]
            if len(c) != 0:
                output_file = ["-o", "".join([x[0] for x in c]) + ".csv"]
            else:
                output_file = ["-o", "none.csv"]

            args = protgraph.parse_args(statistics + configuration + output_file + input_files)
            protgraph.prot_graph(**args)
