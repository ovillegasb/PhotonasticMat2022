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
        "XyzOrdonnee": "1"
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
    lines += "Deltatemps         {:>20.1e}\n".format(options["Deltatemps"] * 1e-12)
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
        self.natypes = 0
        self.atypes = []

        self.dftypes = []

    def write_atominfo(self, MOL, **kwargs):
        lines = ""
        lines_pot = "\n* Fonction potentielle\n"

        for i in MOL.dftypes.index:
            atype = MOL.dftypes.loc[i, "type"]

            if atype not in self.atypes:
                self.atypes.append(atype)
                lines += "\n"
                lines += "* Atome {}\n".format(self.natypes)
                lines += "nom             {:>}\n".format(MOL.dftypes.loc[i, "type"])
                lines += "masse           {:>.2e} kg/mol\n".format(MOL.dftypes.loc[i, "mass"] / 1000)

                lines_pot += "Potentiel {} {}  LJ epsilon {:.7e} J rc {:.1f} - sigma {:.3e} metre\n".format(
                    self.natypes,
                    self.natypes,
                    MOL.dftypes.loc[i, "epsilon"] * 1000 / N_A,  # J
                    self.options["rc"],  # nanometers
                    MOL.dftypes.loc[i, "sigma"] * 1e-10  # metre
                )

                self.natypes += 1
                self.dftypes.append(dict(MOL.dftypes.loc[i, :]))

        self._lines1 += lines
        self._lines_pot += lines_pot

    def write_ffpar(self, **kwargs):
        lines = self._lines0
        lines += "NbTypesAtomes   {:>}\n".format(self.natypes)
        lines += "maille_long     {:.3f} {:.3f} {:.3f} ang\n".format(
            fatomes.options["maille_long"],
            fatomes.options["maille_long"],
            fatomes.options["maille_long"])
        lines += "maille_angle    90.0 90.0 90.0 degre\n"
        lines += "maille_orient   0\n"
        lines += "maille_ref\n"
        lines += "*\n"

        self._lines0 += lines
        self.dftypes = pd.DataFrame(self.dftypes)
        lines_pot = ""
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

        fatomes.lines = self._lines0 + self._lines1 + self._lines_pot

        # writing all
        with open("FAtomes.in", "w") as f:
            f.write(fatomes.lines)

        print("\nfile \033[1;36mFAtomes.in\033[m writed\n")
