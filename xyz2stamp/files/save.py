"""
Submodule dedicated to functions that write information to files.


"""

import time
from scipy.constants import N_A
import itertools as it
import pandas as pd
import numpy as np


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
        "StepLimit": 10,
        "StepAvg": "1",
        "FichierAtomes": "./FAtomes.in",  # config
        "Decoupage": "10 10 10",  # config
        "Parallele": "1 1 1",
        "ConditionX": "Periodique",
        "ConditionY": "Periodique",
        "ConditionZ": "Periodique",
        "drVerlet": "1.e-10",
        "Protection": 10,
        "ProtectionExaStamp": "1",  # config
        "Reprise": "0",  # config
        "TransfertDesForces": "1",
        "RandSeed": "777",
        "Vitesse": "1",
        "NbFantomesSup": "0",
        "XyzSortie": "1",  # config
        "XyzFrequence": 100,
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
    lines += "* Created using XYZ2STAMP module\n"
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

    print("file \033[1;36mDONNEES.in\033[m writed\n")


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

        """
        Se inicializa la clase cargando el sistem y definiendo las opciones.

        """

        self.system = []

        # Head
        now = time.ctime()
        lines0 = "*\n"
        lines0 += "* run file created on {}\n".format(now)
        lines0 += "* Created using XYZ2STAMP module\n"
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

        self.dftypes = []

    def write_atominfo(self, MOL, **kwargs):
        lines = ""
        lines_pot = "* Fonction potentielle\n"
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
                lines += "type            {:>}\n".format("Atome")
                lines += "masse           {:>.2e} kg/mol\n".format(MOL.dftypes.loc[i, "mass"] / 1000)
                if self.natypes == 0:
                    lines += "maille_long     {:.3f} {:.3f} {:.3f} ang\n".format(
                        fatomes.options["maille_long"],
                        fatomes.options["maille_long"],
                        fatomes.options["maille_long"])
                    lines += "maille_angle    90.0 90.0 90.0 degre\n"
                    lines += "maille_orient   0\n"
                    lines += "maille_ref\n"
                    lines += "structure       {:>}\n".format(fatomes.options["structure"])

                lines_pot += "Potentiel {} {}  LJ epsilon {:.7e} J rc {:.1f} - sigma {:.3e} metre\n".format(
                    self.natypes,
                    self.natypes,
                    MOL.dftypes.loc[i, "epsilon"] * 1000 / N_A,  # J
                    self.options["rc"],  # nanometers
                    MOL.dftypes.loc[i, "sigma"] * 1e-10  # metre
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
                        MOL.dfbonds.loc[i, "kb"]
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
                        MOL.dfangles.loc[i, "kth"]
                    )
                    self.angtypes.append(atype)
                    self.n_parintra += 1
        # DIHEDRALS
        if "types" in MOL.dfdih:
            for i in MOL.dfdih.index:
                dtype = MOL.dfdih.loc[i, "types"]
                if dtype not in self.dihtypes and atype[::-1] not in self.dihtypes:
                    lines_ffintra += "{}{:>7}{:>4}{:>4}{:>4} {:>8.2f} kcal/mol {:>8d} {:>8d} {:>8.2f} deg\n".format(
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
                    self.n_parintra += 1

        self._lines1 += lines
        self._lines_pot += lines_pot
        self._lines_xyz += lines_xyz
        self._lines_zmat += lines_zmat
        self._lines_ffintra += lines_ffintra
        self._lines_pcharges += lines_pcharges

    def write_ffpar(self, **kwargs):
        lines = ""
        lines += "NbTypesAtomes   {:>}\n".format(self.natypes)

        self._lines0 += lines
        self.dftypes = pd.DataFrame(self.dftypes)
        lines_pot = ""
        # Add combination using Lorentz-Berthelot rules
        for i, j in it.combinations(range(self.natypes), 2):
            epsilon_ij = np.sqrt(self.dftypes.loc[i, "epsilon"] * self.dftypes.loc[j, "epsilon"])
            sigma_ij = (self.dftypes.loc[i, "sigma"] + self.dftypes.loc[j, "sigma"]) / 2

            lines_pot += "Potentiel {} {}  LJ epsilon {:.7e} J rc {:.1f} - sigma {:.3e} metre\n".format(
                i,
                j,
                epsilon_ij * 1000 / N_A,  # J
                self.options["rc"],  # nanometers
                sigma_ij * 1e-10  # metre
                )

        self._lines_pot += lines_pot

    def write_topol(self, **kwargs):
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
