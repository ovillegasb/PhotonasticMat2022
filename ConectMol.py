#!/bin/env python

"""Script to generate a system molecule from connectivity."""

from stamptools import read_fatomes
from moltools.structure import connectivity
import networkx as nx
import json


fatoms = "FAtomes_000125000.in"
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
    bulk[ipol]["index"] = " ".join(["%d to %d" % (i, i + 9) for i in mol[0:-1:size]])
    ipol += 1
    if len(mol) > atoms:
        atoms = len(mol)
        indexs = mol

# print(indexs)
# print(atoms)

# for vmd
select = " ".join(["%d to %d" % (i, i + 9) for i in indexs[0:-1:10]])
print(select)

#

with open("polysystem.json", "w") as f:
    json.dump(bulk, f, indent=0)

# print(bulk)

"""
init = 0
end = 0

ixat = 0

while ixat < len(indexs):
    at1 = indexs[ixat]
    fxat = 0
    for at2 in indexs[ixat+1:]:
        if at2 - at1 == 1:
            at1 = at2
            fxat = at2
            # print(at1)
            # print(at2)
            # print(fxat)

        elif at2 - at1 > 1:
            ixat = at2 + 1
            ixat = fxat + 1
            print(ixat, fxat)
            break

    # break
"""