#!/bin/env python

"""Module used to generate a gromacs file from a FAtomes file."""

from stamptools import read_fatomes, save_gro
import argparse


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="FAtoms2Gro",
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


def main():
    """Define main function."""
    args = options()

    fatoms = args["fatoms"]
    print(f"FAtomes file: {fatoms}")

    dfatoms, box = read_fatomes(fatoms)
    dfatoms["resid"] = 0
    dfatoms["resname"] = "UKN"

    # This can be changed to be entered as optional parameters.
    atomsRes = 10
    resname = "BUT"

    Natoms = len(dfatoms)
    TotRes = int(Natoms / atomsRes)
    print("TotRes:", TotRes)
    for res in range(TotRes):
        # print("RES:", res)
        for at in range(atomsRes):
            # print("Atom:", at + res*atomsRes)
            dfatoms.loc[at + res*atomsRes, "resid"] = res
            dfatoms.loc[at + res*atomsRes, "resname"] = resname

    save_gro(dfatoms, box)


if __name__ == '__main__':
    main()
