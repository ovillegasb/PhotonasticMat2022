#!/bin/env python

"""
Script used to estimate the time to reach a density using a piston in Stamp.

An input file is required, FAtomes.in, there are two optional arguments, the
DONNEES.in file and the number of steps.

Orlando Villegas
2022.05.11


To receive help use execute:

    python DensPiston.py -h

Examples:

    python DensPiston.py FAtomes.in

    python DensPiston.py FAtomes.in -d DONNEES.in

    python DensPiston.py FAtomes.in -d DONNEES_NPT.in -n 300000

    python DensPiston.py FATOMES/FAtomes_000100000.in -n 200000

"""

from sys import argv
import numpy as np
import pandas as pd
from scipy.constants import N_A
import argparse
import re


def print_mess(message, **kwargs):
    """Print message in color"""
    print("\033[1;36m%s\033[m" % message, **kwargs)

"""Generate command line interface."""
parser = argparse.ArgumentParser(
        prog="DensPiston",
        usage="%(prog)s FAtomes.in [DONNEE.in optional] [-nsteps] value",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"  #, description=__doc__
    )

# FAtomes.in file
parser.add_argument(
    "fatoms",
    help="FAtomes.in file",
    default="",
    type=str
)

# DONNEES.in file
parser.add_argument(
    "-d", "--donnees",
    help="DONNEES.in file, is optional",
    type=str
)

# Nsteps
parser.add_argument(
    "-n", "--nsteps",
    help="Number of steps for piston",
    default=100000,
    type=int
)

args = vars(parser.parse_args())

""" Regular expression that extracts matrix XYZ """
atoms = re.compile(r"""
        ^\s+
        (?P<atsb>[A-Za-z]+\d?\d?)\s+      # Atom name.
        (?P<x>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for X.
        (?P<y>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Y.
        (?P<z>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Z.
        """, re.X)

fatoms = args["fatoms"]
donnees = args["donnees"]
nsteps = args["nsteps"]
decoupage = np.array([1, 1, 1])

print(f"FAtomes file: {fatoms}")
print(f"Nsteps: {nsteps}")

if donnees:
    print(f"DONNEES file: {donnees}")
    with open(donnees) as DAT:
        for line in DAT:
            if "Decoupage" in line:
                decoupage = np.array(line.split()[1:4]).astype(np.int64)
                break

"""READ Fatomes"""
natypes = 0
atomsM = {}
xyz = []
with open(fatoms, "r") as FATM:
    for line in FATM:
        if "NbTypesAtomes" in line:
            line = line.split()
            natypes += int(line[1])
            continue

        if "nom" in line:
            line = line.split()
            nom = line[1]
            continue
            
        if "masse" in line:
            line = line.split()
            mass = np.float64(line[1])
            atomsM[nom] = mass
            continue

        if "maille_long" in line:
            line = line.split()
            box = np.array(line[1:4]).astype(np.float64)
            continue

        if "PositionDesAtomesCart" in line:
            Natoms = int(FATM.readline())
            continue

        if atoms.match(line):
            m = atoms.match(line)
            xyz.append(m.groupdict())
        

print("Number of atoms in XYZ matrix:", Natoms)

print_mess("ATOMS types and Mass [Kg/mol]")
print(atomsM)

box *= decoupage
print_mess("Box dimensions [angs]:")
print(box)

tabXYZ = pd.DataFrame(xyz)

tabXYZ.astype({
    "x": np.float64,
    "y": np.float64,
    "z": np.float64
})

tabXYZ["mass"] = tabXYZ["atsb"].apply(lambda x: atomsM[x])

MassTOT = tabXYZ["mass"].sum() / N_A
MassTOT = MassTOT * decoupage[0] * decoupage[1] * decoupage[2]
print_mess("Total system mass [kg]:", end=" ")
print(MassTOT)

# Box length
Lx = box[0] * 1e-10
Ly = box[1] * 1e-10
Lz = box[2] * 1e-10

x_0 = box[0]

# System volume initial
Vol_0 = Lx * Ly * Lz

print_mess("Volume initial [m^3]:", end=" ")
print(Vol_0)

dens_0 = MassTOT / Vol_0
print_mess("Density initial [kg/m^3]:", end=" ")
print(dens_0)

print("\033[1;33mEnter the desired densitys\033[m")
dens_f = np.float64(input("(in kg/m^3): "))

Vol_f = MassTOT / dens_f
print_mess("Volume final [m^3]:", end=" ")
print(Vol_f)

x_f = Lx + (Vol_f - Vol_0) / (Ly*Lz) # in m
x_f *= 1e10 # to angs
print("The final x-dimension of the box Using the x component [ang]:", x_f)

print("The change of x-component is:", x_f - x_0)

dt = 1 # fs, 1e-15
time = nsteps * dt

# 1e-5 ang / fs
vel = 1e5 * (x_f - x_0) / time
print_mess("The velocity obtained is [1e-5 angs / fs = m / s]: %e" % vel)

print("Divided for each direction (<-->) [m / s]: %e" % (vel / 2))
