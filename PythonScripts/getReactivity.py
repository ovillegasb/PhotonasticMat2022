#!/bin/env python

import argparse
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

def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="getReactivity",
        usage="%(prog)s FAtomes.in [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    # FAtomes.in file
    parser.add_argument(
        "fatoms",
        help="FAtomes.in file",
        default="",
        type=str
    )

    return vars(parser.parse_args())


def read_fatomes(file):
    """Read Fatomes file."""
    natypes = 0
    atomsM = {}
    xyz = []
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

    print("Number of atoms in XYZ matrix:", Natoms)

    print("ATOMS types and Mass [Kg/mol]")
    print(atomsM)

    print("Box dimensions [angs]:")
    print(box)

    tabXYZ = pd.DataFrame(xyz)

    tabXYZ.astype({
        "x": np.float64,
        "y": np.float64,
        "z": np.float64
    })

    tabXYZ["mass"] = tabXYZ["atsb"].apply(lambda x: atomsM[x])

    return tabXYZ, Natoms


def save_reactivity(coords):
    """Save a file reactivity.dat."""

    nat = len(coords)
    name = "reactivity.dat"

    lines = "* ========== \n"
    lines += "* Reactivite \n"
    lines += "* ========== \n"
    lines += "Reactivite \n"
    lines += "%d\n" % nat
    for i in coords.index:
        lines += "%d %s\n" % (i, coords.loc[i, "reactivity"])

    # write lines
    with open(name, "w", encoding="utf-8") as REACT:
        REACT.write(lines)

    print("File reactivity.dat saved")


def main():
    args = options()

    fatoms = args["fatoms"]
    print(f"FAtomes file: {fatoms}")

    dfatoms, Natoms = read_fatomes(fatoms)
    dfatoms["reactivity"] = "INERTE"

    # This can be changed to be entered as optional parameters.
    atomsRes = 10
    i_at_acceptor = 0
    i_at_donor = 3

    for i in range(int(Natoms / atomsRes)):
        dfatoms.loc[i_at_acceptor, "reactivity"] = "ACCEPTEUR"
        dfatoms.loc[i_at_donor, "reactivity"] = "DONNEUR"

        i_at_acceptor += 10
        i_at_donor += 10

    # saving file reactivity
    save_reactivity(dfatoms)


if __name__ == '__main__':
    main()
