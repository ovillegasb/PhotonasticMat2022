#!/usr/bin/env python
# -*- coding=utf-8 -*-

import mdtraj as md
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput
import sys
import os


# Create folder where the files will be saveds
if 'gaussinputs' not in os.listdir():
    os.mkdir('gaussinputs')
else:
    pass

# Molecule resname
res = sys.argv[1]

# Reads the system trajectory
t = md.load('traj_nopbc.xtc', top='confout.gro')

# Every 10 frames
step = 10
traj = t[0:-1:step]

# Extracts the system topology
top = traj.topology
table, bonds = top.to_dataframe()

# The indices of the atoms of interest are selected.
iat_res = top.select("resname %s" % res)

# for frame to gaussian input
for iframe in range(traj.n_frames):
    atoms = table.loc[iat_res, "element"].values
    coord = traj.xyz[:, iat_res][iframe] * 10.0 # to angstroms

    mol = Molecule(atoms, coord)

    # save file gaussian
    GaussianInput(
        mol,
        charge=0,
        spin_multiplicity=1,
        title="Azobenzene Sampling - frame %d" % iframe,
        functional="B3LYP",
        basis_set="Def2SVP",
        route_parameters={
            "SCF": "Tight",
            "TD": "(NStates=3)"
        },
        link0_parameters={
            "%chk": "azoC_%d.chk" % iframe
        }
    ).write_file("gaussinputs/azoC_%d.com" % iframe)