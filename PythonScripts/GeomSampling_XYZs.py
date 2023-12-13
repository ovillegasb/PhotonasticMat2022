#!/usr/bin/env python
# -*- coding=utf-8 -*-

from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput
import sys
import os
import glob
from molcraft.structure import load_xyz


files_xyz = sys.argv()[1]

XYZs = glob.glob("xyz_center_0.3/*.xyz")
print("Number of files:", len(XYZs))

exit()

for conf in XYZs:
    print(conf)
    _, rep, iframe, = conf.split("/")[-1].split("_")
    iframe = iframe.split(".")[0]
    title = "Azobenzene cis Sampling - replica %s frame %s" % (rep, iframe)
    print(title)

    dfatoms = load_xyz(conf)

    atsb = list(dfatoms.loc[0:25, "atsb"].values)
    xyz = dfatoms.loc[0:25, ["x", "y", "z"]].values

    mol = Molecule(atsb, xyz)

    # save file gaussian
    GaussianInput(
        mol,
        charge=0,
        spin_multiplicity=1,
        title=title,
        functional="Cam-B3LYP",
        basis_set="6-311+g(d,p)",
        route_parameters={
            "TD": "(NStates=6)"
        },
        link0_parameters={
            "%chk": "azoC_%s_%s.chk" % (rep, iframe)
        }
    ).write_file("azoC_%s_%s.com" % (rep, iframe))

    break


"""
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
"""
