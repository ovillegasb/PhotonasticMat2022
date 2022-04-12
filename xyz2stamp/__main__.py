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
    parser.add_argument("-t", "--type",
                        help="donnees or fatomes",
                        default="donnees",
                        metavar="entrytype",
                        type=str,
                        choices=["donnees", "fatomes"])

    donnees = parser.add_argument_group("donnees file parameters")

    # donnees group
    donnees.add_argument("-ensemble", "--Ensemble",
                        help="Thermodynamic Ensemble",
                        default="NVT",
                        type=str,
                        choices=["NVT", "NPT"])


    return vars(parser.parse_args())


def main():
    args = options()

    save.write_run_file(**args)



# RUN

main()