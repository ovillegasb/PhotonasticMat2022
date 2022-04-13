"""
Module created to generate input files to STAMP.

Author: Orlando VILLEGAS
Date: 2022
\033[1;36m
 __   ____     ___________   _____ _______       __  __ _____
 \\ \\ / /\\ \\   / /___  /__ \\ / ____|__   __|/\\   |  \\/  |  __ \\
  \\ V /  \\ \\_/ /   / /   ) | (___    | |  /  \\  | \\  / | |__) |
   > <    \\   /   / /   / / \\___ \\   | | / /\\ \\ | |\\/| |  ___/
  / . \\    | |   / /__ / /_ ____) |  | |/ ____ \\| |  | | |
 /_/ \\_\\   |_|  /_____|____|_____/   |_/_/    \\_\\_|  |_|_|
\033[m
Module created to generate input files to STAMP.

Author: Orlando VILLEGAS
Date: 2022

###############################################################

"""

import argparse
from xyz2stamp.files import save
from xyz2stamp.files.load import load_structure
from xyz2stamp.structure import connectivity, MOL, FField


TITLE = """
Module created to generate input files to STAMP.

Author: Orlando VILLEGAS
Date: 2022
\033[1;36m
 __   ____     ___________   _____ _______       __  __ _____
 \\ \\ / /\\ \\   / /___  /__ \\ / ____|__   __|/\\   |  \\/  |  __ \\
  \\ V /  \\ \\_/ /   / /   ) | (___    | |  /  \\  | \\  / | |__) |
   > <    \\   /   / /   / / \\___ \\   | | / /\\ \\ | |\\/| |  ___/
  / . \\    | |   / /__ / /_ ____) |  | |/ ____ \\| |  | | |
 /_/ \\_\\   |_|  /_____|____|_____/   |_/_/    \\_\\_|  |_|_|
\033[m
Module created to generate input files to STAMP.

Author: Orlando VILLEGAS
Date: 2022

###############################################################

"""


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


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

    donnees = parser.add_argument_group(
        "\033[1;36mDonnees file parameters\033[m")

    # donnees group
    donnees.add_argument(
        "-ensemble", "--Ensemble",
        help="Thermodynamic Ensemble",
        default="NVT",
        type=str,
        choices=["NVT", "NPT"]
    )

    # Temperature
    donnees.add_argument(
        "-temp", "--Trequis",
        help="Temperature (K)",
        default=300.0,
        type=float,
        metavar="TEMPERATURE"
    )

    # Time step
    donnees.add_argument(
        "-dt", "--Deltatemps",
        help="Time step (fs)",
        default=1.0,
        type=float,
        metavar="TIME"
    )

    # Number of steps
    donnees.add_argument(
        "-nstep", "--StepLimit",
        help="Numer of steps N x fs = Simulation time.",
        default=10,
        type=int,
        metavar="NSTEP"
    )

    # save frame every step
    donnees.add_argument(
        "-frames", "--XyzFrequence",
        help="save a frame every n step.",
        default=10,
        type=int,
        metavar="NSTEPS"
    )

    # save a checkpoint every n step
    donnees.add_argument(
        "-check", "--Protection",
        help="save a checkpoint every n step.",
        default=10,
        type=int,
        metavar="NSTEPS"
    )

    fatomes = parser.add_argument_group(
        "\033[1;36mFatomes file parameters\033[m")

    # fatomes group
    fatomes.add_argument(
        "-f", "--files",
        help="Specifies all coordinates files in series that you want to load",
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

    print(TITLE)

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
            print("No file found")

        # save.write_topol(**args)


# RUN

main()
