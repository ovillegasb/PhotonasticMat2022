#!/bin/env python

"""Script to generate a system molecule from connectivity."""

from stamptools.analysis import read_fatomes
from moltools.structure import connectivity
import networkx as nx
import json
import argparse


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="ConectMol",
        usage="%(prog)s FAtomes.in [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    # FAtomes.in file
    parser.add_argument(
        "fatoms",
        help="FAtomes.in file",
        default="",
        type=str
    )

    return vars(parser.parse_args())


def group_runs(li, tolerance=2):
    """Group consecutive integers and tolerate gaps of 1."""
    out = []
    last = li[0]
    for x in li:
        if x-last > tolerance:
            yield out
            out = []
        out.append(x)
        last = x
    yield out


args = options()

fatoms = args["fatoms"]
print(f"FAtomes file: {fatoms}")

dfatoms, box, connects = read_fatomes(fatoms)

print(dfatoms)
# print(connects)

conn = connectivity()

conn.define_atoms(dfatoms)

conn.read_dict(connects)

# print(conn)
# print(conn[0])
# print(conn.subgraph([0]))
# print(nx.node_connected_component(conn, 0))

atoms_MOL = nx.weakly_connected_components(conn)

# First molecule
# print(next(atoms_MOL))

size = 10  # esto se tiene que cambiar

# searching most big polymer
atoms = 0
indexs = []

# Dict n atoms : indexs
bulk = dict()
ipol = 0
for mol in atoms_MOL:
    mol = list(sorted(mol))
    bulk[ipol] = dict()
    bulk[ipol]["Natoms"] = len(mol)
    # bulk[ipol]["index"] = " ".join(["%d to %d" % (i, i + 9) for i in mol[0:-1:size]])
    # print(sorted(mol))
    m = sorted(mol)
    groups = group_runs(m)

    # print(list(groups))

    ndx = " ".join(["%d to %d" % (lm[0], lm[-1]) for lm in groups])
    bulk[ipol]["index"] = ndx
    # bulk[ipol]["index"] = " ".join([str(i) for i in mol])
    ipol += 1
    if len(mol) > atoms:
        atoms = len(mol)
        indexs = mol


with open("polysystem.json", "w") as f:
    json.dump(bulk, f, indent=0, sort_keys=True)

print("\nFile polysystem.json saved.")
