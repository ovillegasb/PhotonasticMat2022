#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""Generates a GRO trajectory of the centers of mass of the molecules."""

import sys
import os
import pandas as pd
import numpy as np
import mdtraj as md
from stamptools.stamp import STAMP
from molcraft.structure import Elements


try:
    sysType = sys.argv[1]
    donnees = sys.argv[2]
except IndexError:
    print("Error")
    print("CMtraj [STAMP or GRO] [donnees/gro] [optional traj] [optional resid\
from 0]")
    exit()

try:
    resid_select = int(sys.argv[4])
except IndexError:
    resid_select = None

if sysType == "GRO":
    top = sys.argv[2]
    traj = sys.argv[3]


def save_gro(coord, box, title="GRO FILE", time=0.0, out="."):
    """
    Save an gro file of coordinates.

    Parameters:
    -----------
    coord : DataFrame

    name : str

    box = array(1x3)

    """
    nat = len(coord)
    lines = ""
    lines += "{}, t= {:.3f}\n".format(title, time)
    lines += "%5d\n" % nat

    for i in coord.index:
        lines += "{:>8}{:>7}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
            str(coord.loc[i, "resid"]) + coord.loc[i, "resname"],
            coord.loc[i, "atsb"],
            i + 1,
            coord.loc[i, "x"] * 0.1,  # + box[0] / 20
            coord.loc[i, "y"] * 0.1,  # + box[1] / 20
            coord.loc[i, "z"] * 0.1)  # + box[2] / 20

    lines += "   {:.5f}   {:.5f}   {:.5f}\n".format(
        box[0] / 10,
        box[1] / 10,
        box[2] / 10)

    return lines


def center_of_mass(coords, box, masses):
    """Compute the center of mass, the mass weighterd barycenter."""
    theta_i = coords / box * 2 * np.pi
    xi_i = np.cos(theta_i)
    eta_i = np.sin(theta_i)
    xi_m = np.sum(xi_i * masses[:, np.newaxis], axis=0) / masses.sum()
    eta_m = np.sum(eta_i * masses[:, np.newaxis], axis=0) / masses.sum()
    theta_m = np.arctan2(-eta_m, -xi_m) + np.pi

    return box * theta_m / 2 / np.pi


output = "traj_cm.gro"
conf = "confout_cm.gro"
lines = ""
begin = 0

if os.path.exists(output):
    print(f"Trajectory file {output} exist")
    traj_cm = md.load(output, top=conf)
    print("Number of frames analysed:", traj_cm.n_frames)
    begin = traj_cm.n_frames - 1

else:
    with open(output, "w") as GRO:
        GRO.write(lines)

if sysType == "STAMP":
    system = STAMP(donnees)
    time_per_frame = system.time_per_frame
    box_in_frame = system.box_in_frame
    atoms_per_mol = system.atoms_per_mol
    connectivity = system.connectivity
    top = system.topology
    # Load trajectory
    traj = system.get_traj()
    resid = 0
    #####
    dt = np.float64(system.donnees['Deltatemps']) * int(system.donnees["XyzFrequence"]) * 1e12
    time = 0.0
    #####
    for n_frame, frame in enumerate(traj):
        if n_frame < begin:
            continue
        traj_cm = []
        box = box_in_frame[n_frame][0:3]
        ## time = time_per_frame.loc[n_frame, "time"]
        time = n_frame * dt
        for imol, mol in enumerate(atoms_per_mol):
            if resid_select is not None:
                if imol != resid_select:
                    continue
            mol_ndx = atoms_per_mol[mol]
            masses = top.loc[mol_ndx["index"], "mass"].values
            coords = frame.loc[mol_ndx["index"], :]
            resname = "POL"
            symbol = "C"
            if imol == 0:
                resname = "PHO"
                symbol = "P"

            mol_cm = center_of_mass(
                coords.loc[:, ["x", "y", "z"]].values,
                box,
                masses
            )
            traj_cm.append({
                "atsb": symbol,
                "resname": resname,
                "resid": imol + 1,
                "x": mol_cm[0],
                "y": mol_cm[1],
                "z": mol_cm[2]
                })

        traj_cm = pd.DataFrame(traj_cm)
        lines += save_gro(traj_cm, box, time=time)
        if n_frame == 0:
            with open(conf, "w") as GRO:
                GRO.write(lines)

        with open(output, "a") as GRO:
            GRO.write(lines)
        lines = ""

elif sysType == "GRO":
    # Reads the system trajectory
    traj = md.load(traj, top=top)

    # Extracts the system topology
    top = traj.topology
    table, bonds = top.to_dataframe()
    boxs = traj.unitcell_lengths
    dt = traj.timestep
    time = 0.0

    nmol = len(table["resSeq"].unique())
    for n_frame, frame in enumerate(traj):
        if n_frame < begin:
            continue
        print("Frame:", n_frame, "Total:", len(traj))
        traj_cm = []
        box = boxs[n_frame] * 10
        time = n_frame * dt
        for mol in range(nmol):
            if resid_select is not None:
                if mol != resid_select:
                    continue

            # The indices of the atoms of interest are selected.
            iat_res = top.select("resid %s" % mol)
            atoms = table.loc[iat_res, "element"]
            masses = atoms.apply(lambda x: Elements[x]["mass"]).values
            coords = traj.xyz[:, iat_res][n_frame] * 10.0  # to angstroms
            resname = "THF"
            symbol = "O"
            if mol == 0:
                resname = "PHO"
                symbol = "P"
            
            mol_cm = center_of_mass(
                coords,
                box,
                masses
            )

            traj_cm.append({
                "atsb": symbol,
                "resname": resname,
                "resid": mol + 1,
                "x": mol_cm[0],
                "y": mol_cm[1],
                "z": mol_cm[2]
                })

        traj_cm = pd.DataFrame(traj_cm)
        lines += save_gro(traj_cm, box, time=time)
        if n_frame == 0:
            with open(conf, "w") as GRO:
                GRO.write(lines)
        with open(output, "a") as GRO:
            GRO.write(lines)

        lines = ""

else:
    print("Error")
    print("CMtraj [STAMP or GRO] [donnees/gro] [optional traj] [optional resid\
from 0]")
    exit()
