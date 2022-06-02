
"""
Module that stores different functions and tools for working with STAMP files.

Orlando Villegas - 2022

"""

import re
import numpy as np
import pandas as pd

""" Regular expression that extracts matrix XYZ """
atoms = re.compile(r"""
        ^\s+
        (?P<atsb>[A-Za-z]+\d?\d?)\s+      # Atom name.
        (?P<x>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for X.
        (?P<y>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Y.
        (?P<z>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Z.
        """, re.X)


def read_fatomes(file):
    """Read Fatomes file."""
    natypes = 0
    atomsM = {}
    xyz = []
    connects = dict()
    with open(file, "r") as FATM:
        for line in FATM:
            if "*" == line[0]:
                # ignore lines with the * symbol
                continue

            elif "NbTypesAtomes" in line:
                line = line.split()
                natypes += int(line[1])
                continue

            elif "nom" in line:
                line = line.split()
                nom = line[1]
                continue

            elif "masse" in line:
                line = line.split()
                mass = np.float64(line[1])
                atomsM[nom] = mass
                continue

            elif "maille_long" in line:
                line = line.split()
                box = np.array(line[1:4]).astype(np.float64)
                continue

            elif "PositionDesAtomesCart" in line:
                Natoms = int(FATM.readline())
                continue

            elif atoms.match(line):
                m = atoms.match(line)
                xyz.append(m.groupdict())

            elif "Zmatrice" in line:
                N = int(FATM.readline())
                print("N conectivity:", N)
                for _ in range(N):
                    zline = FATM.readline()
                    zline = zline.split()
                    zline = [int(i) for i in zline]
                    connects[zline[0]] = zline[1:]

    print("Number of atoms in XYZ matrix:", Natoms)

    print("ATOMS types and Mass [Kg/mol]")
    print(atomsM)

    print("Box dimensions [angs]:")
    print(box)

    tabXYZ = pd.DataFrame(xyz)

    tabXYZ = tabXYZ.astype({
        "x": np.float64,
        "y": np.float64,
        "z": np.float64
    })

    tabXYZ["mass"] = tabXYZ["atsb"].apply(lambda x: atomsM[x])

    return tabXYZ, box, connects


def save_gro(table, box, name="coord.gro"):
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
        box[0] * 0.1,
        box[1] * 0.1,
        box[2] * 0.1))

    GRO.close()
    print("\nSaved gro file: \033[1;36m%s\033[m writed\n" % gro)
