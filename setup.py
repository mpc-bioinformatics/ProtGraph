from setuptools import find_packages, setup

# Get Packages from requirements.txt by parsing them
packages = []
with open("requirements.txt", "r") as reqs:
    packages = reqs.readlines()

# read README.md
with open("README.md", "r", encoding="utf-8") as long_desc:
    long_description = long_desc.read()

setup(
    name='protgraph',
    version='0.3.11',
    author="Dominik Lux",
    description="ProtGraph, a graph generator for proteins.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mpc-bioinformatics/ProtGraph",
    project_urls={
        "Bugs": "https://github.com/mpc-bioinformatics/ProtGraph/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: BSD License',
    ],
    license="BSD",
    python_requires=">=3.7",
    entry_points=dict(console_scripts=[
        'protgraph = protgraph.protgraph:main',
        'protgraph_pepsqlite_to_fasta = protgraph.scripts.pepsqlite_to_fasta:main [sqlite]',
        'protgraph_replace_fasta_header = protgraph.scripts.replace_fasta_header:main',
        'protgraph_generate_fasta_decoys = protgraph.scripts.generate_fasta_decoys:main',
        'protgraph_compact_fasta = protgraph.scripts.compact_fasta:main',
        'protgraph_print_sums = protgraph.scripts.print_sums:main',
        'protgraph_convert_fasta_to_sp_embl_txt = protgraph.scripts.convert_fasta_to_sp_embl_txt:main'
    ]),
    packages=find_packages(),
    include_package_data=True,
    install_requires=packages,
    extras_require={
        "postgres": ["psycopg>=3.0"],
        "mysql": ["mysql"],
        "sqlite": ["apsw>=3.42.0.0"],
        "cassandra": ["cassandra-driver"],
        "gremlin": ["gremlinpython"],
        "redis": ["redis", "redisgraph"],
        "full": ["mysql", "psycopg>=3.0", "apsw>=3.42.0.0", "cassandra-driver", "redis", "redisgraph", "gremlinpython"],
    },
)
