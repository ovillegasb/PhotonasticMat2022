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
import time
import itertools as it
import pandas as pd
import numpy as np
from moltools.structure import MOL, ATOM
from moltools.ffield import get_atoms_types, get_ffparameters, get_interactions_list
from mkitpgmx import save_gro, save_itp

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

References:


Docs:

http://ambermd.org/vdwequation.pdf

"""


def write_run(**kwargs):
    """
    Reads the options set to generate the donnees file.

    """

    # set the default values
    options = {
        "Ensemble": "NVT",
        "Trequis": 100.0,
        "TauNVT": "3.e-12",
        "Schema": "verlet",
        "Deltatemps": 1.0,
        "StepLimit": 20,
        "StepAvg": "1",
        "FichierAtomes": "./FAtomes.in",  # config
        "Decoupage": "5 5 5",  # config
        "Parallele": "2 2 2",
        "ConditionX": "Periodique",
        "ConditionY": "Periodique",
        "ConditionZ": "Periodique",
        "drVerlet": "2.e-10",
        "Protection": 10,
        "ProtectionExaStamp": "1",  # config
        "Reprise": "0",  # config
        "TransfertDesForces": "1",
        "RandSeed": "777",
        "Vitesse": "1",
        "NbFantomesSup": "0",
        "XyzSortie": "1",  # config
        "XyzFrequence": 10,
        "XyzOrdonnee": "1",
        "Molecule": 1,
        "Molecule_CalculIntra": 1,
        "Molecule_ContribVdwIntra12": 0,
        "Molecule_ContribVdwIntra13": 0,
        "Molecule_ContribVdwIntra14": 0
    }

    # update the default values based on the arguments
    options.update(kwargs)

    # Head
    now = time.ctime()
    lines = "*\n"
    lines += "* run file created on {}\n".format(now)
    lines += "* Created using MOL2STAMP module\n"
    lines += "* mail: orlando.villegas@chimieparistech.psl.eu\n"
    lines += "*\n"
    lines += "Ensemble           {:>20}\n".format(options["Ensemble"])
    lines += "Trequis            {:>20.1f}\n".format(options["Trequis"])
    lines += "TauNVT             {:>20}\n".format(options["TauNVT"])
    lines += "*\n"
    lines += "Schema             {:>20}\n".format(options["Schema"])
    lines += "Deltatemps         {:>20.1e}\n".format(options["Deltatemps"] * 1e-15)
    lines += "StepLimit          {:>20d}\n".format(options["StepLimit"])
    lines += "StepAvg            {:>20}\n".format(options["StepAvg"])
    lines += "*\n"
    lines += "FichierAtomes      {:>20}\n".format(options["FichierAtomes"])
    lines += "Decoupage          {:>20}\n".format(options["Decoupage"])
    lines += "Parallele          {:>20}\n".format(options["Parallele"])
    lines += "ConditionX         {:>20}\n".format(options["ConditionX"])
    lines += "ConditionY         {:>20}\n".format(options["ConditionY"])
    lines += "ConditionZ         {:>20}\n".format(options["ConditionZ"])
    lines += "*\n"
    lines += "drVerlet           {:>20}\n".format(options["drVerlet"])
    lines += "*\n"
    lines += "Protection         {:>20d}\n".format(options["Protection"])
    lines += "ProtectionExaStamp {:>20}\n".format(options["ProtectionExaStamp"])
    lines += "Reprise            {:>20}\n".format(options["Reprise"])
    lines += "*\n"
    lines += "TransfertDesForces {:>20}\n".format(options["TransfertDesForces"])
    lines += "*\n"
    lines += "RandSeed           {:>20}\n".format(options["RandSeed"])
    lines += "Vitesse            {:>20}\n".format(options["Vitesse"])
    lines += "NbFantomesSup      {:>20}\n".format(options["NbFantomesSup"])
    lines += "*\n"
    lines += "XyzSortie          {:>20}\n".format(options["XyzSortie"])
    lines += "XyzFrequence       {:>20d}\n".format(options["XyzFrequence"])
    lines += "XyzOrdonnee        {:>20}\n".format(options["XyzOrdonnee"])
    lines += "*\n"
    lines += "Molecule                     {:>20d}\n".format(options["Molecule"])
    lines += "Molecule_CalculIntra         {:>20d}\n".format(options["Molecule_CalculIntra"])
    lines += "Molecule_ContribVdwIntra12   {:>20d}\n".format(options["Molecule_ContribVdwIntra12"])
    lines += "Molecule_ContribVdwIntra13   {:>20d}\n".format(options["Molecule_ContribVdwIntra13"])
    lines += "Molecule_ContribVdwIntra14   {:>20d}\n".format(options["Molecule_ContribVdwIntra14"])

    # writing all
    with open("DONNEES.in", "w") as f:
        f.write(lines)

    print("\nfile \033[1;36mDONNEES.in\033[m writed\n")


class fatomes:

    options = {
        "NbTypesAtomes": 0,
        "structure": "FICHIER",
        "maille_long": 12.4,  # angs
        "maille_angle": 90.0,
        "maille_orient": 0,
        "maille_red": None,
        "rc": 1.5
    }

    lines = ""

    def __init__(self):

        """The class is initialized by loading the system and defining the options."""

        self.system = []

        # Head
        now = time.ctime()
        lines0 = "*\n"
        lines0 += "* run file created on {}\n".format(now)
        lines0 += "* Created using MOL2STAMP module\n"
        lines0 += "* mail: orlando.villegas@chimieparistech.psl.eu\n"
        lines0 += "*\n"

        self._lines0 = lines0
        self._lines1 = ""
        self._lines_pot = ""
        self._lines_xyz = "*\nPositionDesAtomesCart angstrom\nNATOMS\n"
        self._lines_zmat = "*\n* Connectivity\nZmatrice\nNATOMS\n"
        self._lines_ffintra = "*\n* Champ de force intramoleculaire\nChampDeForces\nNPARINTRA\n"
        self._lines_pcharges = "*\nModificationChargeDesAtomes e-\nNATOMS\n"

        self.natypes = 0
        self.natoms = 0
        self.n_parintra = 0
        self.atypes = []
        self.btypes = []
        self.angtypes = []
        self.dihtypes = []
        self.imptypes = []
        self.dftypes = []

    def write_atominfo(self, MOL, **kwargs):
        lines = ""
        lines_pot = "*\n* Fonction potentielle\n"
        lines_xyz = ""
        lines_zmat = ""
        lines_ffintra = ""
        lines_pcharges = ""
        # ATOMS TYPES
        for i in MOL.dftypes.index:
            atype = MOL.dftypes.loc[i, "type"]

            if atype not in self.atypes:
                self.atypes.append(atype)
                lines += ""
                lines += "* Atome {}\n".format(self.natypes)
                lines += "nom             {:>}\n".format(MOL.dftypes.loc[i, "type"])
                lines += "nomFF           {:>}\n".format(MOL.dftypes.loc[i, "type"])
                lines += "nomXYZ          {:>}\n".format(MOL.dftypes.loc[i, "type"])
                lines += "type            {:>}\n".format("Atome")
                lines += "masse           {:>.2e} kg/mol\n".format(MOL.dftypes.loc[i, "mass"] / 1000)
                if self.natypes == 0:
                    lines += "structure       {:>}\n".format(fatomes.options["structure"])
                    lines += "maille_long     {:.3f} {:.3f} {:.3f} ang\n".format(
                        fatomes.options["maille_long"],
                        fatomes.options["maille_long"],
                        fatomes.options["maille_long"])
                    lines += "maille_angle    90.0 90.0 90.0 degre\n"
                    lines += "maille_orient   0\n"
                    lines += "maille_ref\n"

                lines_pot += "Potentiel {} {}  LJ epsilon {:.4f} kcal/mol rc {:.1f} - sigma {:.4f} ang\n".format(
                    self.natypes,
                    self.natypes,
                    MOL.dftypes.loc[i, "epsilon"],       # kcal/mol
                    self.options["rc"],                  # nanometers
                    MOL.dftypes.loc[i, "sigma"]          # angstrom
                )

                self.natypes += 1
                self.dftypes.append(dict(MOL.dftypes.loc[i, :]))

            self.natoms += 1

            lines_xyz += "{:>6}{:>15.6f}{:>15.6f}{:>15.6f}\n".format(
                MOL.dftypes.loc[i, "type"],
                MOL.dftypes.loc[i, "x"],
                MOL.dftypes.loc[i, "y"],
                MOL.dftypes.loc[i, "z"]
            )

            neighbors = list(MOL.connect.neighbors(i))
            neighbors = [str(i) for i in neighbors]

            lines_zmat += "{} {}\n".format(i, ' '.join(neighbors))
            lines_pcharges += "{}  {:>8.3f}\n".format(i, MOL.dftypes.loc[i, "charge"])
        # BONDS
        if "types" in MOL.dfbonds:
            for i in MOL.dfbonds.index:
                btypes = MOL.dfbonds.loc[i, "types"]
                if btypes not in self.btypes:
                    lines_ffintra += "{}{:>8}{:>4} {:>8.3f} ang {:>8.2f} kcal/mol/ang2\n".format(
                        "bond_gaff",
                        btypes[0],
                        btypes[1],
                        MOL.dfbonds.loc[i, "b0"],
                        MOL.dfbonds.loc[i, "kb"] * 2
                    )
                    self.btypes.append(btypes)
                    self.n_parintra += 1
        # ANGLES
        if "types" in MOL.dfangles:
            for i in MOL.dfangles.index:
                atype = MOL.dfangles.loc[i, "types"]
                if atype not in self.angtypes and atype[::-1] not in self.angtypes:
                    lines_ffintra += "{}{:>7}{:>4}{:>4} {:>8.2f} deg {:>8.2f} kcal/mol/rad2\n".format(
                        "angle_gaff",
                        atype[0],
                        atype[1],
                        atype[2],
                        MOL.dfangles.loc[i, "th0"],
                        MOL.dfangles.loc[i, "kth"] * 2
                    )
                    self.angtypes.append(atype)
                    self.n_parintra += 1
        # DIHEDRALS
        if "types" in MOL.dfdih:
            for i in MOL.dfdih.index:
                dtype = MOL.dfdih.loc[i, "types"]
                if dtype not in self.dihtypes and atype[::-1] not in self.dihtypes:
                    lines_ffintra += "{}{:>7}{:>4}{:>4}{:>4} {:>8.2f} kcal/mol {:>8d} - {:>8d} - {:>8.2f} deg\n".format(
                        "torsion_gaff",
                        dtype[0],
                        dtype[1],
                        dtype[2],
                        dtype[3],
                        MOL.dfdih.loc[i, "Vn"],
                        MOL.dfdih.loc[i, "divider"],
                        MOL.dfdih.loc[i, "n"],
                        MOL.dfdih.loc[i, "phi"]
                    )
                    self.dihtypes.append(dtype)
                    self.dihtypes.append(dtype[::-1])
                    self.n_parintra += 1
        # IMPROPERS
        if "types" in MOL.dfimp:
            for i in MOL.dfimp.index:
                itype = MOL.dfimp.loc[i, "types"]
                if itype not in self.imptypes:
                    lines_ffintra += "{}{:>7}{:>4}{:>4}{:>4} {:>8.2f} kcal/mol\n".format(
                        "impropre_gaff",
                        itype[2],  # For stamp, cetral atom is the first
                        itype[0],
                        itype[1],
                        itype[3],
                        MOL.dfimp.loc[i, "Vn"]
                    )
                    self.imptypes.append(itype)
                    self.n_parintra += 1

        self._lines1 += lines
        self._lines_pot += lines_pot
        self._lines_xyz += lines_xyz
        self._lines_zmat += lines_zmat
        self._lines_ffintra += lines_ffintra
        self._lines_pcharges += lines_pcharges

    def _write_ffpar(self, **kwargs):
        lines = ""
        lines += "NbTypesAtomes   {:>}\n".format(self.natypes)

        self._lines0 += lines
        self.dftypes = pd.DataFrame(self.dftypes)
        lines_pot = ""
        # Add combination using Lorentz-Berthelot rules
        for i, j in it.combinations(range(self.natypes), 2):
            epsilon_ij = np.sqrt(self.dftypes.loc[i, "epsilon"] * self.dftypes.loc[j, "epsilon"])
            sigma_ij = (self.dftypes.loc[i, "sigma"] + self.dftypes.loc[j, "sigma"]) / 2

            lines_pot += "Potentiel {} {}  LJ epsilon {:.4f} kcal/mol rc {:.1f} - sigma {:.4f} ang\n".format(
                i,
                j,
                epsilon_ij,          # kcal/mol
                self.options["rc"],  # nanometers
                sigma_ij             # ang
                )

        self._lines_pot += lines_pot

    def write_topol(self, **kwargs):
        self._write_ffpar(**kwargs)
        fatomes.lines += self._lines0
        fatomes.lines += self._lines1
        fatomes.lines += self._lines_pot
        self._lines_xyz = self._lines_xyz.replace('NATOMS', '%d' % self.natoms)
        fatomes.lines += self._lines_xyz
        self._lines_pcharges = self._lines_pcharges.replace('NATOMS', '%d' % self.natoms)
        fatomes.lines += self._lines_pcharges
        self._lines_zmat = self._lines_zmat.replace('NATOMS', '%d' % self.natoms)
        fatomes.lines += self._lines_zmat
        self._lines_ffintra = self._lines_ffintra.replace('NPARINTRA', '%d' % self.n_parintra)
        fatomes.lines += self._lines_ffintra

        # writing all
        with open("FAtomes.in", "w") as f:
            f.write(fatomes.lines)

        print("\nfile \033[1;36mFAtomes.in\033[m writed\n")


def options():
    """Generate command line interface."""

    parser = argparse.ArgumentParser(
        prog="XYZ2STAMP",
        usage="%(prog)s [-option] value",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"  #, description=__doc__
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
        choices=["NVE", "NVT", "NPT"]
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


def print_steps(message):
    print("\033[1;35m%s\033[m" % message)


def main():
    """
    Central core of program execution.

    """

    print(TITLE)

    args = options()

    if args["type"] == "donnees":
        # Generates the DONNEES input file with the simulation parameters.
        write_run(**args)

    elif args["type"] == "fatomes":
        # Generates the input file FAtoms with the structural and force
        # field information of the system.

        if args["files"]:
            # 0) The molecule class is initialized.
            mol = MOL()

            # 1) Read all entry files.
            print_steps("0) Read all entry files.")
            print("files:", args["files"])
            print("N files", len(args["files"]))

            # INIT for ...
            mol.load_file(args["files"][0])
            print("DATA:\n", mol.dfatoms)

            # 2) Search connectivity from geometry.
            print_steps("2) Search connectivity from geometry.")
            mol.search_connectivity()

            # 3) Atom types are assigned.
            print_steps("3) Atom types are assigned.")
            get_atoms_types(mol, args["forcefield"])

            # 4) The force field parameters are assigned.
            print_steps("4) The force field parameters are assigned.")
            get_interactions_list(mol)
            get_ffparameters(mol, args["forcefield"])

            Fatomes = fatomes()
            Fatomes.write_atominfo(mol)
            # END for ...

            # 5) write all information.
            print_steps("5) write all information.")
            Fatomes.write_topol()
        else:
            print("No file found")


if __name__ == "__main__":
    # RUN

    main()
