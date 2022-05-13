
"""
Script used to estimate the time to reach a density using a piston in Stamp.

Two input files FAtomes and PasDeCalcul_Iteration_XXXXXXX.xyz are required.

Orlando Villegas
2022.05.11

Steps
-----

1. The dimensions of the box are read from the xyz file. This version assumes a
rectangular box.

2. The volume of the box is calculated.

3. Read the types of atoms, the masses.

"""

from sys import argv
import numpy as np
import pandas as pd
from scipy.constants import N_A

xyz = argv[1]
fatoms = argv[2]

print(f"XYZ file: {xyz}")
print(f"FAtomes file: {fatoms}")

with open(xyz, "r") as XYZ:
    # for line in XYZ:

    # First line
    # number of atoms
    Natoms = int(XYZ.readline())

    box = np.array(XYZ.readline().split()[:3]).astype(np.float64)

natypes = 0
atomsM = {}
with open(fatoms, "r") as FATM:
    # lines = FATM.readlines()
    for line in FATM:
        if "NbTypesAtomes" in line:
            natypes += int(line.split()[1])
            continue
            
        if natypes != 0:
            line = line.split()
            if "nom" in line:
                nom = line[1]
                atomsM
            elif "masse" in line:
                mass = np.float64(line[1])
                atomsM[nom] = mass

tabXYZ = pd.read_csv(
    xyz,
    sep=r"\s+",
    header=None,
    skiprows=2,
    names=["atsb", "x", "y", "z", "0", "1"]
)

tabXYZ.astype({
    "x": np.float64,
    "y": np.float64,
    "z": np.float64
})


print("Number od atoms in file XYZ:", Natoms)

print("Box dimensions")
print(box)

print("ATOMS types and Mass [Kg/mol]")
print(atomsM)

tabXYZ["mass"] = tabXYZ["atsb"].apply(lambda x: atomsM[x])
MassTOT = tabXYZ["mass"].sum() / N_A
print("Total system mass [kg]:", MassTOT)

Lx = box[0] * 1e-10
Ly = box[1] * 1e-10
Lz = box[2] * 1e-10

x_0 = box[0]

Vol_0 = Lx * Ly * Lz
print("Volume initial [m^3]:", Vol_0)

dens_0 = MassTOT / Vol_0
print("Density initial [kg/m^3]:", dens_0)

print("Enter the desired density")
dens_f = np.float64(input("(in kg/m^3): "))

Vol_f = MassTOT / dens_f
print("Volume final [m^3]:", Vol_f)

x_f = Lx + (Vol_f - Vol_0) / (Ly*Lz) # in m
x_f *= 1e10 # to angs
print("The final x-dimension of the box Using the x component [ang]:", x_f)

print("The change of x-component is:", x_f - x_0)
print(r"5 % of change:", (x_f - x_0) * 5 / 100)

dt = 1 # fs, 1e-15
steps= 5000 # for 5ps for steps from piston
time = steps * dt
# 1e-5 ang / fs
vel = 1e5 * (x_f - x_0) * 5 / 100 / time
print("The velocity obtained is [1e-5 angs / fs = m / s]: %e" % vel)

print("Divided for each direction (<-->) [m / s]: %e" % (vel / 2))
print("The total number of steps required (for each step of 1fs) is: %d" % int((x_f - x_0)/(1e-5*vel)))
