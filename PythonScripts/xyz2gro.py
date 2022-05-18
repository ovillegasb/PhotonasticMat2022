#!/bin/env python

"""

Script used to convert STAMP xyz files to GROMACS gro files.

Inputs options:
    -r/--res : [N1-RES1-atomsRES1, N2-RES2-atomsRES2, ...]
    -p/--vpiston : velocity piston

"""

import argparse
import pandas as pd
import numpy as np


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
            prog="xyz2gro",
            usage="%(prog)s files.xyz -res N1-RES1-atomsRES1 ... [-p] velocity positive(m/s)",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="Enjoy the program!"  #, description=__doc__
        )

    """Add the arguments."""

    # file.xyz file
    parser.add_argument(
        "xyz",
        help="file.xyz",
        default="",
        type=str
    )

    # resname options
    parser.add_argument(
        "-r", "--res",
        help="Resname options, N1-RES1-atomsRES1 N2-RES2-atomsRES2 ...",
        type=str,
        nargs="+"
    )

    # piston option
    parser.add_argument(
        "-p", "--vpiston",
        help="Piston velocity (m/s)",
        type=float,
        default=None
    )

    return vars(parser.parse_args())


def save_gro(table, name, box):
    """Save coordinate to file *.gro from dataframe with x, y, z."""
    nat = len(table)
    gro = name

    GRO = open(gro, "w", encoding="utf-8")
    GRO.write("GRO FILE\n")
    GRO.write("%5d\n" % nat)
    for i in table.index:
        GRO.write("{:>8}{:>7}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
            str(table.loc[i, "resid"]) + table.loc[i, "resname"],
            table.loc[i, "atsb"],
            i + 1,
            table.loc[i, "x"] * 0.1,
            table.loc[i, "y"] * 0.1,
            table.loc[i, "z"] * 0.1)
        )
    GRO.write("   {:.5f}   {:.5f}   {:.5f}\n".format(
        box[0] / 10,
        box[1] / 10,
        box[2] / 10))

    GRO.close()
    print("\nSaved gro file: \033[1;36m%s\033[m writed\n" % gro)

args = options()


xyz = args["xyz"]
nameXYZ = xyz.split(".")[-2]
step = int(nameXYZ.split("_")[-1])
res = args["res"]
gro = nameXYZ + ".gro"

print("File xyz:", xyz)
print("Step:", step)

print("Number od residues:", len(res))
print("Output file:", gro)

velp = 0.0
if args["vpiston"]:
    print("Piston is activated")
    velp += args["vpiston"]*1e-5
    print("Piston velocity [angs/fs]:", velp)

Natoms = 0
with open(xyz, "r") as XYZ:
    # for line in XYZ:
    #     print(line)

    # first line is number od atoms
    line = XYZ.readline()
    Natoms += int(line)

    # second line, box information
    # assuming a rectagular box
    line = XYZ.readline()
    box = np.array(line.split()[0:3]).astype(np.float64)

print("Number of atoms in file:", Natoms)
print("BOx information (angstroms):", box)

coord = pd.read_csv(
    xyz,
    sep=r"\s+",
    header=None,
    names=["atsb", "x", "y", "z", "o1", "o2"],
    skiprows=2
)

coord.astype({
    "x": np.float64,
    "y": np.float64,
    "z": np.float64
})


# coordGRO = coord.copy()

"""
natomstot = 0
for rmol in res:
    nmolr, mol, natomsr = rmol.split("-")
    nmolr = int(nmolr)
    natomsr = int(natomsr)
    natomstot += natomsr*nmolr
    print(nmolr)
    print(mol)
    print(natomsr)
""" 


Lx, Ly, Lz = box
natomstot = 0
atoms_added = 0
ires = 1
ires_added = 0
RES = []
idRES = []
for i in range(len(res)):
    rmol = res[i]
    nmolr, mol, natomsr = rmol.split("-")
    nmolr = int(nmolr)
    natomsr = int(natomsr)
    natomstot += natomsr*nmolr
    # print("RES:", mol)
    # print("NMOL:", nmolr)
    # print("NATOMS:", natomsr)
    # print("NATOMS TOTAL:", natomstot)

    while ires <= nmolr + ires_added:
        for ati in range(atoms_added, natomsr + atoms_added):
            coord.loc[ati, "x"] = coord.loc[ati, "x"] + Lx/2 + (step*velp)/2
            coord.loc[ati, "y"] = coord.loc[ati, "y"] + Ly/2
            coord.loc[ati, "z"] = coord.loc[ati, "z"] + Lz/2
            idRES.append(ires)
            RES.append(mol)
            atoms_added += 1
        ires += 1
    ires_added += nmolr

coord["resname"] = RES
coord["resid"] = idRES

save_gro(coord, name=gro, box=box)
