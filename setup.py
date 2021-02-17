from setuptools import find_packages, setup

# Get Packages from Pipfile by parsing them
lines = []
with open("Pipfile", "r") as pipfile:
    lines = pipfile.readlines()
c = [idx for idx, x in enumerate(lines) if x[0] == "["]

# dev packages index range
dev_packages_sec = [x for x in c if "[dev-packages]" in lines[x]][0]
if c.index(dev_packages_sec) + 1 == len(c):
    dev_packages_end = len(lines)
else:
    dev_packages_end = c[c.index(dev_packages_sec)+1]

# packages index range
packages_sec = [x for x in c if "[packages]" in lines[x]][0]
if c.index(packages_sec) + 1 == len(c):
    packages_end = len(lines)
else:
    packages_end = c[c.index(packages_sec)+1]

# retrieve the packages (except itself)
packages = []
for i in range(dev_packages_sec+1, dev_packages_end):
    if lines[i] == "\n":
        continue
    package = lines[i].split("=")[0].strip()
    if package != "protgraph":
        packages.append(package)
for i in range(packages_sec+1, packages_end):
    if lines[i] == "\n":
        continue
    package = lines[i].split("=")[0].strip()
    if package != "protgraph":
        packages.append(package)


setup(
    name='protgraph',
    version='0.0.1a',
    description="ProtGraph, a graph generator for proteins.",
    python_requires=">=3.6",
    entry_points=dict(console_scripts=['protgraph=protgraph:main']),
    packages=find_packages(),
    install_requires=packages
)
