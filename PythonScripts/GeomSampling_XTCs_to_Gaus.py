#!/usr/bin/env python
# -*- coding=utf-8 -*-

import numpy as np
import mdtraj as md
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput
import sys
import os

out = "SamplingUV-Vis"

try:
    top = sys.argv[1]
    trj = sys.argv[2]
    b = int(sys.argv[3])
    e = int(sys.argv[4])
    out += "_{}_{}".format(sys.argv[3], sys.argv[4])
except IndexError:
    print("Options: .gro .xtc begin end")
    exit()

if not os.path.exists(out):
    os.mkdir(out)

resid = 0  # start from 1
replica = "thf"
isomer = input("Isomer, cis or trans: ")

print("File top:", top)
print("File trj:", trj)
print("Time >=", b, "ps")
print("Time <", e, "ps")
print("Resid:", resid)

# Reads the system trajectory
t = md.load(trj, top=top)

# Every 10 frames
step = None
traj = t[b:e]

frames_sample = np.random.choice(range(len(traj)), 100, replace=False)

# Extracts the system topology
top = traj.topology
table, bonds = top.to_dataframe()

# The indices of the atoms of interest are selected.
iat_res = top.select("resid %s" % resid)

# for frame to gaussian input
for i in frames_sample:
    atoms = table.loc[iat_res, "element"].values
    coord = traj.xyz[:, iat_res][i] * 10.0  # to angstroms

    mol = Molecule(atoms, coord)

    # save file gaussian
    GaussianInput(
        mol,
        charge=0,
        spin_multiplicity=1,
        title="Azobenzene %s Sampling - sample %d" % (isomer, i),
        functional="Cam-B3LYP",
        basis_set="6-311+g(d,p)",
        route_parameters={
            "TD": "(NStates=6)",
            "SCRF": "(Solvent=TetraHydroFuran)"
        },
        link0_parameters={
            "%mem": "8GB",
            "%nprocshared": "12"
        }
    ).write_file("%s/azo%s_%s_%005d.com" % (out, isomer[0].upper(), replica, i))
