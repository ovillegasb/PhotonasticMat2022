"""
Submodule dedicated to functions that write information to files.


"""

import time
from scipy.constants import N_A


def write_run(**kwargs):
    """
    Reads the options set to generate the donnees file.

    """

    # set the default values
    options = {
        "Ensemble": "NVT",  # config
        "Trequis": "100.",  # config
        "TauNVT": "3.e-12",
        "Schema": "verlet",
        "Deltatemps": "1.0e-15",  # config
        "StepLimit": "10",  # config
        "StepAvg": "1",
        "FichierAtomes": "./FAtomes.in",  # config
        "Decoupage": "10 10 10",  # config
        "Parallele": "1 1 1",
        "ConditionX": "Periodique",
        "ConditionY": "Periodique",
        "ConditionZ": "Periodique",
        "drVerlet": "1.e-10",
        "Protection": "10",  # config
        "ProtectionExaStamp": "1",  # config
        "Reprise": "0",  # config
        "TransfertDesForces": "1",
        "RandSeed": "777",
        "Vitesse": "1",
        "NbFantomesSup": "0",
        "XyzSortie": "1",  # config
        "XyzFrequence": "100",  # config
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
    lines += "Trequis            {:>20}\n".format(options["Trequis"])
    lines += "TauNVT             {:>20}\n".format(options["TauNVT"])
    lines += "*\n"
    lines += "Schema             {:>20}\n".format(options["Schema"])
    lines += "Deltatemps         {:>20}\n".format(options["Deltatemps"])
    lines += "StepLimit          {:>20}\n".format(options["StepLimit"])
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
    lines += "Protection         {:>20}\n".format(options["Protection"])
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
    lines += "XyzFrequence       {:>20}\n".format(options["XyzFrequence"])
    lines += "XyzOrdonnee        {:>20}\n".format(options["XyzOrdonnee"])

    # writing all
    with open("DONNEES.in", "w") as f:
        f.write(lines)

    print("file DONNEES.in writed")


def fatomes(MOL):

    def write_topol(**kwargs):
        # Head
        now = time.ctime()
        lines = "*\n"
        lines += "* run file created on {}\n".format(now)
        lines += "* Created using XYZ2STAMP module\n"
        lines += "* mail: orlando.villegas@chimieparistech.psl.eu\n"
        lines += "*\n"
        lines += "NbTypesAtomes      {:>20}\n".format(kwargs["NbTypesAtomes"])

        for i in MOL.dftypes.index:
            lines += "* Atome {}\n".format(i)
            lines += "nom             {}\n".format(MOL.dftypes.loc[i, "type"])
            lines += "type            {}\n".format(kwargs["type"])
            lines += "masse           {:.2e} kg/mol\n".format(MOL.dftypes.loc[i, "mass"] / 1000)
            lines += "structure       {}\n".format(kwargs["structure"])
            lines += "* Box\n"
            lines += "maille_long     {:.3e} {:.3e} {:.3e} metre\n".format(kwargs["box_a"], kwargs["box_a"], kwargs["box_a"])
            lines += "maille_angle    90.0 90.0 90.0 degre\n"
            lines += "maille_orient   0\n"
            lines += "maille_ref\n"
            lines += "*\n"
            lines += "* Fonction potentielle\n"
            lines += "Potentiel {} {}  LJ epsilon {:.7e} J rc {:.1f} - sigma {:.3e} metre\n".format(
                    i,
                    i,
                    MOL.dftypes.loc[i, "epsilon"] * 1000 / N_A,  # J
                    kwargs["rc"],  # nanometers
                    MOL.dftypes.loc[i, "sigma"] * 1e-10  # metre
                )
        """
        Potentiel 0 0  LJ epsilon 1.6567944e-21 J rc 2.5 - sigma 3.405e-10 metre
        """
        # writing all
        with open("FAtomes.in", "w") as f:
            f.write(lines)

    options = {
        "NbTypesAtomes": len(MOL.dftypes),
        "type": "Atome",
        "structure": "CFC",
        "box_a": 3.0 * 1e-10,  # Angstrom
        "rc": 2.0  # nanometers
    }

    write_topol(**options)

    print("file FAtomes.in writed")
