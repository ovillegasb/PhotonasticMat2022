#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""Sampling."""

import sys
import os
import numpy as np
from stamptools.stamp import STAMP
from stamptools.analysis import center_of_mass
# from molcraft.structure import save_xyz
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

out = "SamplingUV-Vis"

try:
    file = sys.argv[1]
except IndexError:
    print("Dont file DONNEES.in")
    exit()

# b = input("begin time in ps: "
b = float(input("begin time in ps: "))
e = float(input("end time in ps: "))
resid = 0
replica = int(input("Replica number: "))
isomer = input("Isomer, cis or trans: ")

print("File:", file)
print("Time >=", b, "ps")
print("Time <", e, "ps")
print("Resid:", resid)
out += "_{}_{}".format(str(int(b)), str(int(e)))
print("out folder:", out)
if not os.path.exists(out):
    os.mkdir(out)

system = STAMP(file)

atoms_per_mol = system.atoms_per_mol
connectivity = system.connectivity
top = system.topology
box = system.box
# vol = system.vol
time_per_frame = system.time_per_frame
b = time_per_frame[time_per_frame["time"] >= b].index[0]
e = time_per_frame[time_per_frame["time"] < e].index[-1]

mol_ndx = atoms_per_mol[resid]

print("Initial Frame:", b)
print("Final Frame :", e)

traj = system.get_traj(b=b, e=e)

frames_sample = np.random.choice(range(len(traj)), 100, replace=False)
# 2print(len(frames_sample))

for i in frames_sample:
    xyz = traj[i]
    name = "mol_%d_%005d" % (replica, i)
    print(name)
    title = "Azobenzene %s Sampling - replica %s sample %s" % (isomer, replica, i)
    mol_xyz = xyz.loc[mol_ndx["index"], :]
    mol_conn = connectivity.sub_connect(mol_ndx["index"])
    # update coordinates
    mol_conn.update_coordinates(mol_xyz)
    # remove PBC
    mol_conn.noPBC(box, center=np.zeros(3))
    # Reset index and symbols, and add mass
    mol_conn = mol_conn.reset_nodes()
    mol_conn.simple_at_symbols(add_mass=True)
    # Search and add hydrogen to vacant atoms
    mol_conn.add_hydrogen(box, type_add="terminal")
    new_mol_xyz = mol_conn.get_df()

    cm = center_of_mass(
        new_mol_xyz.loc[:, ["x", "y", "z"]].values,
        new_mol_xyz.loc[:, "mass"].values
    )

    new_mol_xyz["x"] -= cm[0]
    new_mol_xyz["y"] -= cm[1]
    new_mol_xyz["z"] -= cm[2]

    atsb = list(new_mol_xyz.loc[:, "atsb"].values)
    xyz = new_mol_xyz.loc[:, ["x", "y", "z"]].values
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
        link0_parameters={  # "%chk": "azoC_%s_%s.chk" % (rep, iframe),
            "%mem": "8GB",
            "%nprocshared": "12"

        }
    ).write_file("%s/azo%s_%d_%005d.com" % (out, isomer[0].upper(), replica, i))

    # save_xyz(new_mol_xyz, name=f"{out}/{name}")
    # exit()
