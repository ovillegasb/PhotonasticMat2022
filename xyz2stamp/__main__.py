r"""
Module created to generate input files to STAMP.

Author: Orlando VILLEGAS
Date: 2022

 __   ____     ___________   _____ _______       __  __ _____
 \ \ / /\ \   / /___  /__ \ / ____|__   __|/\   |  \/  |  __ \
  \ V /  \ \_/ /   / /   ) | (___    | |  /  \  | \  / | |__) |
   > <    \   /   / /   / / \___ \   | | / /\ \ | |\/| |  ___/
  / . \    | |   / /__ / /_ ____) |  | |/ ____ \| |  | | |
 /_/ \_\   |_|  /_____|____|_____/   |_/_/    \_\_|  |_|_|

Module created to generate input files to STAMP.

Author: Orlando VILLEGAS
Date: 2022

###############################################################

"""

import argparse
from xyz2stamp.files import save
from xyz2stamp.files.load import load_structure
from xyz2stamp.structure import connectivity, MOL, FField


def options():
    """Generate command line interface."""

    parser = argparse.ArgumentParser(
        prog="XYZ2STAMP",
        usage="%(prog)s [-option] value",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!",
        description=__doc__
    )

    # Add the arguments

    # Inputfile type
    parser.add_argument(
        "-t", "--type",
        help="donnees or fatomes",
        default="donnees",
        metavar="entrytype",
        type=str,
        choices=["donnees", "fatomes"]
    )

    donnees = parser.add_argument_group("donnees file parameters")

    # donnees group
    donnees.add_argument(
        "-ensemble", "--Ensemble",
        help="Thermodynamic Ensemble",
        default="NVT",
        type=str,
        choices=["NVT", "NPT"]
    )

    fatomes = parser.add_argument_group("fatomes file parameters")

    # fatomes group
    fatomes.add_argument(
        "-f", "--files",
        help="Indica all coordinates files seguidas\
        que quiere cargar",
        default=None,
        nargs="+"
    )

    # Force file
    fatomes.add_argument(
        "-ff", "--forcefield",
        help="Force field family",
        default="gaff",
        type=str,
        choices=["gaff", "oplsaa", "gromos"]
    )

    return vars(parser.parse_args())


def main():
    """
    Central core of program execution.

    """

    args = options()

    if args["type"] == "donnees":
        # Generates the DONNEES input file with the simulation parameters.
        save.write_run(**args)

    elif args["type"] == "fatomes":
        # Generates the input file FAtoms with the structural and force
        # field information of the system.

        if args["files"]:
            print(args["files"])
            print(len(args["files"]))
            # for ...
            mol0 = args["files"][0]

            dfatoms = load_structure(mol0)
            print(dfatoms)

            connect = connectivity()
            connect.get_connectivity(dfatoms)

            df = connect.get_df()
            print(df)

            print("Natoms", connect.number_of_nodes())
            print("Nbonds", connect.number_of_edges())

            print("Edges", connect.edges)

            MOL(dfatoms, connect)

            # print(connect[0])
            # print(connect[1])

            FField(args["forcefield"])

            print(MOL.dftypes)

            save.fatomes(MOL)

        else:
            print("anyone")

        # save.write_topol(**args)


# RUN

main()
