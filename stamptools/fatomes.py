"""Module to define the FATOME class."""

from stamptools.analysis import decoTime
from stamptools.analysis import center_of_mass_polar, translate_to
from molcraft import structure
import itertools as it
import numpy as np
import pandas as pd
from datetime import datetime
import re
import os


# dat location
location = os.path.dirname(os.path.realpath(__file__))
cell_unit = os.path.join(location, "oplsaa.dat")


# Some definitions

NAMES_OPLS = [
    "CT", "CM", "HT", "HM", "OHopls", "HOopls"
]

NAMES_PC = [
    "ca", "ha", "oh", "ho", "ne", "nf", "cb"
]

ATOMSTYPES_par = {
    "DU": {
        "nomXYZ": "O",
        "nomFF": "DU",
        "type": "Atome",
        "masse": "1.0000000000e-06 kg/mol",
        "charge": "0.0000000000e+00 e-",
        "gele": "1"
    },
    "CT": {
        "nomXYZ": "C",
        "nomFF": "CT",
        "type": "Atome",
        "masse": "1.2011000000e-02 kg/mol",
        "charge": "-1.2000000000e-01 e-"
    },
    "CM": {
        "nomXYZ": "C",
        "nomFF": "CM",
        "type": "Atome",
        "masse": "1.2011000000e-02 kg/mol",
        "charge": "-1.1500000000e-01 e-"
    },
    "HT": {
        "nomXYZ": "H",
        "nomFF": "HT",
        "type": "Atome",
        "masse": "1.0008000000e-03 kg/mol",
        "charge": "6.0000000000e-02 e-"
    },
    "HM": {
        "nomXYZ": "H",
        "nomFF": "HM",
        "type": "Atome",
        "masse": "1.0008000000e-03 kg/mol",
        "charge": "1.1500000000e-01 e-"
    },
    "ca": {
        "nomXYZ": "C",
        "nomFF": "ca",
        "type": "Atome",
        "masse": "12.e-3 kg/mol"
    },
    "ha": {
        "nomXYZ": "H",
        "nomFF": "ha",
        "type": "Atome",
        "masse": "1.e-3  kg/mol"
    },
    "oh": {
        "nomXYZ": "O",
        "nomFF": "oh",
        "type": "Atome",
        "masse": "16.e-3 kg/mol"
    },
    "ho": {
        "nomXYZ": "H",
        "nomFF": "ho",
        "type": "Atome",
        "masse": "1.e-3  kg/mol"
    },
    "cb": {
        "nomXYZ": "C",
        "nomFF": "ca",
        "type": "Atome",
        "masse": "12.e-3 kg/mol"
    },
    "ne": {
        "nomXYZ": "N",
        "nomFF": "ne",
        "type": "Atome",
        "masse": "14.e-3 kg/mol"
    },
    "nf": {
        "nomXYZ": "N",
        "nomFF": "nf",
        "type": "Atome",
        "masse": "14.e-3 kg/mol"
    }
}

POTENTIAL_par_OPLS = {
    "OHopls": "LJ sigma 3.12 ang epsilon 0.066 kcal/mol rc 12.5 ang",
    "HOopls": "LJ sigma 0.00 ang epsilon 0.000 kcal/mol rc 12.5 ang",
    "DU": "LJ sigma   3.3996695084e+00 ang epsilon   3.7318093892e-03 eV rc   2.5000000000e+00 -",
    "CT": "LJ sigma 3.50 ang epsilon 0.066 kcal/mol rc 12.5 ang",
    "CM": "LJ sigma 3.55 ang epsilon 0.076 kcal/mol rc 12.5 ang",
    "HT": "LJ sigma 2.50 ang epsilon 0.030 kcal/mol rc 12.5 ang",
    "HM": "LJ sigma 2.42 ang epsilon 0.030 kcal/mol rc 12.5 ang"
}

POTENTIAL_par_GAFF = {
    "ca": "LJ sigma   3.3996695084e+00 ang epsilon   3.7318093892e-03 eV rc   2.5000000000e+00 -",
    "ha": "LJ sigma   2.5996424595e+00 ang epsilon   6.5089698650e-04 eV rc   2.5000000000e+00 -",
    "oh": "LJ sigma   3.0664733878e+00 ang epsilon   9.1299150639e-03 eV rc   2.5000000000e+00 -",
    "ho": "LJ sigma   0.0000000000e+00 ang epsilon   0.0000000000e+00 eV rc   2.5 ang",
    "ne": "LJ sigma   3.2499985238e+00 ang epsilon   7.3768325136e-03 eV rc   2.5000000000e+00 -",
    "nf": "LJ sigma   3.2499985238e+00 ang epsilon   7.3768325136e-03 eV rc   2.5000000000e+00 -",
    "cb": "LJ sigma   3.3996695084e+00 ang epsilon   3.7318093892e-03 eV rc   2.5000000000e+00 -"
}

BONDS_par_OPLS = {
    ("HOopls", "OHopls"): "bond_opls     HOopls  OHopls     0.945 ang   553.0 kcal/mol/ang2   0.0 -  0.0 -",
    ("CT", "OHopls"): "bond_opls     CT  OHopls     1.410 ang   320.0 kcal/mol/ang2   0.0 -  0.0 -",
    ("CT",  "CT"): "bond_opls     CT  CT     1.529 ang       268.0 kcal/mol/ang2   0.0 -   0.0 -",
    ("CT",  "HT"): "bond_opls     CT  HT     1.090 ang       340.0 kcal/mol/ang2   0.0 -   0.0 -",
    ("CT",  "CM"): "bond_opls     CT  CM     1.510 ang       317.0 kcal/mol/ang2   0.0 -   0.0 -",
    ("CM",  "HM"): "bond_opls     CM  HM     1.080 ang       340.0 kcal/mol/ang2   0.0 -   0.0 -",
    ("CM",  "CM"): "bond_opls     CM  CM     1.340 ang       549.0 kcal/mol/ang2   0.0 -   0.0 -"
}

BONDS_par_GAFF = {
    ("ca", "ca"): "bond_gaff ca ca       1.3870000000 ang     478.4000000000 kcal/mol/ang2",
    ("ca", "ha"): "bond_gaff ca ha       1.0870000000 ang     344.3000000000 kcal/mol/ang2",
    ("ca", "oh"): "bond_gaff ca oh       1.3620000000 ang     386.1000000000 kcal/mol/ang2",
    ("ca", "ne"): "bond_gaff ca ne       1.4310000000 ang     361.8000000000 kcal/mol/ang2",
    ("ca", "nf"): "bond_gaff ca nf       1.4310000000 ang     361.8000000000 kcal/mol/ang2",
    ("oh", "ho"): "bond_gaff oh ho       0.9740000000 ang     369.6000000000 kcal/mol/ang2",
    ("ne", "nf"): "bond_gaff ne nf       1.2475 ang     797.8847 kcal/mol/ang2"
}

ANGLES_par_OPLS = {
    ("CT", "OHopls", "HOopls"): "angle_opls     CT  OHopls  HOopls     108.5 degre  55.0  kcal/mol   0.0 -  0.0 -",
    ("CM", "CT", "OHopls"): "angle_opls     CM  CT  OHopls     109.5 degre  50.0  kcal/mol   0.0 -  0.0 -",
    ("HT", "CT", "OHopls"): "angle_opls     HT  CT  OHopls     109.5 degre  35.0  kcal/mol   0.0 -  0.0 -",
    ("HT", "CT", "HT"): "angle_opls     HT  CT  HT     107.8 degre   33.0  kcal/mol          0.0 -   0.0 -",
    ("CT", "CT", "HT"): "angle_opls     CT  CT  HT     110.7 degre   37.5  kcal/mol          0.0 -   0.0 -",
    ("HT", "CT", "CM"): "angle_opls     HT  CT  CM     109.5 degre   35.0  kcal/mol          0.0 -   0.0 -",
    ("CT", "CM", "HM"): "angle_opls     CT  CM  HM     117.0 degre   35.0  kcal/mol          0.0 -   0.0 -",
    ("CT", "CM", "CM"): "angle_opls     CT  CM  CM     124.0 degre   70.0  kcal/mol          0.0 -   0.0 -",
    ("HM", "CM", "CM"): "angle_opls     HM  CM  CM     120.0 degre   35.0  kcal/mol          0.0 -   0.0 -",
    ("CT", "CT", "CM"): "angle_opls     CT  CT  CM     111.1 degre   63.0  kcal/mol          0.0 -   0.0 -"
}

ANGLES_par_GAFF = {
    ("ca", "ca", "ca"): "angle_gaff ca ca ca     119.9700000000 deg      67.1800000000 kcal/mol",
    ("ca", "ca", "ha"): "angle_gaff ca ca ha     120.0100000000 deg      48.4600000000 kcal/mol",
    ("ca", "ca", "oh"): "angle_gaff ca ca oh     119.9400000000 deg      69.8500000000 kcal/mol",
    ("ca", "ca", "ne"): "angle_gaff ca ca ne     103.4 deg      53.2182 kcal/mol",
    ("ca", "ca", "nf"): "angle_gaff ca ca nf     103.4 deg      53.2182 kcal/mol",
    ("ca", "oh", "ho"): "angle_gaff ca oh ho     109.4700000000 deg      48.8500000000 kcal/mol",
    ("ca", "ne", "nf"): "angle_gaff ca ne nf     110.1 deg      76.4997 kcal/mol",
    ("ca", "nf", "ne"): "angle_gaff ca nf ne     110.1 deg      76.4997 kcal/mol"
}

# https://pubs-acs-org.inc.bib.cnrs.fr/doi/suppl/10.1021/ja9621760/suppl_file/ja11225.pdf
DIHEDRALS_par_OPLS = {
    ("CM", "CT", "OHopls", "HOopls"): "torsion_opls     CM  CT  OHopls  HOopls  -0.356 kcal/mol -0.174  kcal/mol   0.492  kcal/mol",  # amberTools 051
    ("HT", "CT", "OHopls", "HOopls"): "torsion_opls     HT  CT  OHopls  HOopls   0.000 kcal/mol  0.000  kcal/mol   0.450  kcal/mol",  # oplsaa.dat
    ("CM", "CM", "CT", "OHopls"): "torsion_opls     CM  CM  CT  OHopls   1.711 kcal/mol -0.500  kcal/mol   0.663  kcal/mol",
    ("HM", "CM", "CT", "OHopls"): "torsion_opls     HM  CM  CT  OHopls   0.000 kcal/mol  0.000  kcal/mol   0.468  kcal/mol",
    ("HT",  "CT",  "CT",  "HT"): "torsion_opls    HT  CT  CT  HT      0.000  kcal/mol 0.000  kcal/mol         0.300  kcal/mol",
    ("HT",  "CT",  "CT",  "CM"): "torsion_opls    HT  CT  CT  CM      0.000  kcal/mol 0.000  kcal/mol         0.366  kcal/mol",
    ("CT",  "CT",  "CM",  "HM"): "torsion_opls    CT  CT  CM  HM      0.000  kcal/mol 0.000  kcal/mol         0.000  kcal/mol",
    ("HT",  "CT",  "CM",  "HM"): "torsion_opls    HT  CT  CM  HM      0.000  kcal/mol 0.000  kcal/mol         0.318  kcal/mol",
    ("CM",  "CM",  "CT",  "CT"): "torsion_opls    CM  CM  CT  CT      0.346  kcal/mol 0.405  kcal/mol         -0.940  kcal/mol",
    ("CT",  "CM",  "CM",  "CT"): "torsion_opls    CT  CM  CM  CT      0.000  kcal/mol 1400.00  kcal/mol         0.000  kcal/mol",
    ("HM",  "CM",  "CM",  "HM"): "torsion_opls    HM  CM  CM  HM      0.000  kcal/mol 1400.00  kcal/mol         0.000  kcal/mol",
    ("CT",  "CM",  "CM",  "HM"): "torsion_opls    CT  CM  CM  HM      0.000  kcal/mol 1400.00  kcal/mol         0.000  kcal/mol",
    ("HT",  "CT",  "CM",  "CM"): "torsion_opls    HT  CT  CM  CM      0.000  kcal/mol 0.000  kcal/mol         -0.372  kcal/mol",
    ("CM",  "CT",  "CT",  "CM"): "torsion_opls    CM  CT  CT  CM      1.300  kcal/mol -0.05  kcal/mol         0.200  kcal/mol"
}

DIHEDRALS_par_GAFF = {
    ("ca", "ca", "ca", "ca"): "torsion_gaff ca ca ca ca      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ca", "ca", "ha"): "torsion_gaff ca ca ca ha      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ca", "ca", "oh"): "torsion_gaff ca ca ca oh      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ca", "ca", "ne"): "torsion_gaff ca ca ca ne      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ca", "ca", "nf"): "torsion_gaff ca ca ca nf      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ca", "oh", "ho"): "torsion_gaff ca ca oh ho       1.8000000000 kcal/mol   2.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ne", "nf", "ca", "1"): "torsion_gaff ca ne nf ca       17.6364 kcal/mol   1.00 -  2.00 -     180.0 deg",
    ("ca", "ne", "nf", "ca", "2"): "torsion_gaff ca ne nf ca        4.7014 kcal/mol   1.00 -  1.00 -       0.0 deg",
    ("ha", "ca", "ca", "ha"): "torsion_gaff ha ca ca ha      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ha", "ca", "ca", "oh"): "torsion_gaff ha ca ca oh      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ha", "ca", "ca", "ne"): "torsion_gaff ha ca ca ne      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ha", "ca", "ca", "nf"): "torsion_gaff ha ca ca nf      14.5000000000 kcal/mol   4.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ca", "nf", "ne"): "torsion_gaff ca ca nf ne       2.39005736   kcal/mol   1.00 -   2.00 -     180.0000000000 deg",
    ("ca", "ca", "ne", "nf"): "torsion_gaff ca ca ne nf       2.39005736   kcal/mol   1.00 -   2.00 -     180.0000000000 deg"
}

IMPROPERS_par_GAFF = {
    ("ca", "ca", "ca", "ha"): "impropre_gaff    ca ca ca ha          1.1000000000 kcal/mol",
    ("ca", "ca", "ca", "oh"): "impropre_gaff    ca ca ca oh          1.1000000000 kcal/mol",
    ("ca", "ca", "ca", "ne"): "impropre_gaff    ca ca ca ne          1.1000000000 kcal/mol",
    ("ca", "ca", "ca", "nf"): "impropre_gaff    ca ca ca nf          1.1000000000 kcal/mol"
}

ContribDispRepIntra_PC = """28
cb cb 0.0 0.0 0.0 0.0 0.0 0.0
ca ca 0.0 0.0 0.50 0.0 0.0 0.8333
ca ha 0.0 0.0 0.50 0.0 0.0 0.8333
ca oh 0.0 0.0 0.50 0.0 0.0 0.8333
ca ho 0.0 0.0 0.50 0.0 0.0 0.8333
ca ne 0.0 0.0 0.50 0.0 0.0 0.8333
ca nf 0.0 0.0 0.50 0.0 0.0 0.8333
ca cb 0.0 0.0 0.50 0.0 0.0 0.8333
ha ha 0.0 0.0 0.50 0.0 0.0 0.8333
ha oh 0.0 0.0 0.50 0.0 0.0 0.8333
ha ho 0.0 0.0 0.50 0.0 0.0 0.8333
ha ne 0.0 0.0 0.50 0.0 0.0 0.8333
ha nf 0.0 0.0 0.50 0.0 0.0 0.8333
ha cb 0.0 0.0 0.50 0.0 0.0 0.8333
oh oh 0.0 0.0 0.50 0.0 0.0 0.8333
oh ho 0.0 0.0 0.50 0.0 0.0 0.8333
oh ne 0.0 0.0 0.50 0.0 0.0 0.8333
oh nf 0.0 0.0 0.50 0.0 0.0 0.8333
oh cb 0.0 0.0 0.50 0.0 0.0 0.8333
ho ho 0.0 0.0 0.50 0.0 0.0 0.8333
ho ne 0.0 0.0 0.50 0.0 0.0 0.8333
ho nf 0.0 0.0 0.50 0.0 0.0 0.8333
ho cb 0.0 0.0 0.50 0.0 0.0 0.8333
ne ne 0.0 0.0 0.50 0.0 0.0 0.8333
ne nf 0.0 0.0 0.50 0.0 0.0 0.8333
ne cb 0.0 0.0 0.50 0.0 0.0 0.8333
nf nf 0.0 0.0 0.50 0.0 0.0 0.8333
nf cb 0.0 0.0 0.50 0.0 0.0 0.8333"""

Keywords = {}


def change_atsb(x):
    """Change atom types to simple symbols."""
    if x in ["ca", "cb", "CT", "CM", "C", "CTO"]:
        return "C"
    elif x in ["ha", "ho", "HT", "HM", "H", "HOopls", "H1O"]:
        return "H"
    elif x in ["nf", "ne", "N"]:
        return "N"
    elif x in ["oh", "O", "OHopls"]:
        return "O"
    elif x in ["DU"]:
        return "DU"
    else:
        print(f"WARNING: {x} atom not reconize.")
        return x


""" Regular expression for FAtomes """
noms = re.compile(r"""
    nom\s+(?P<nom>\w+)
    """, re.X)

masses = re.compile(r"""
    masse\s+(?P<mass>\d+\.\w+[+-]?\d+)
    """, re.X)

charges = re.compile(r"""
    charge\s+(?P<charge>[+-]?\d+\.\d+\w?[+-]?\d*)
    """, re.X)

""" Regular expression that extracts matrix XYZ """
atoms = re.compile(r"""
        ^\s*
        (?P<atypes>\w{1,6})\s+            # Atom name.
        (?P<x>[+-]?\d+\.\d+\w?[+-]?\d*)\s+      # Orthogonal coordinates for X.
        (?P<y>[+-]?\d+\.\d+\w?[+-]?\d*)\s+      # Orthogonal coordinates for Y.
        (?P<z>[+-]?\d+\.\d+\w?[+-]?\d*)\s+      # Orthogonal coordinates for Z.
        """, re.X)


def combRules_geom(par1, par2):
    """Calculate the VDW parameters using the geometric combination rule."""
    return np.sqrt(par1[0] * par2[0]), np.sqrt(par1[1] * par2[1])


class FATOME:
    """Object representing the STAMP topology file."""

    def __init__(self, name="FAtomes.in"):
        """Initialize the class."""
        self.name = name

    def read_connect(self, BulkConnect):
        """
        Read The connectivity of a system.

        Parameters
        ----------
        BulkConnect : molcraft.structure.BULK
        System connectivity.

        """
        self.bulk = BulkConnect

    def get_atomstypes(self):
        """Generate a dictionary of atom types."""
        atoms_types = {}

        coords = self.bulk.dfatoms.copy()
        box = self.bulk.box
        # print("{} {} {}".format(*tuple(box * 10)))
        # print("{:.6f} {:.6f} {:.6f} ang".format(*tuple(box * 10)))
        print(coords)
        atypes = coords["name"].unique()
        print("Atoms types:", " ".join(atypes))

        for i, atp in enumerate(atypes):
            atoms_types[i] = {}
            atoms_types[i]["nom"] = atp
            for key in ATOMSTYPES_par[atp]:
                atoms_types[i][key] = ATOMSTYPES_par[atp][key]

            if i == 0:
                atoms_types[i]["structure"] = "FICHIER"
                atoms_types[i]["maille_long"] = "{:.6f} {:.6f} {:.6f} ang".format(*tuple(box * 10))
                atoms_types[i]["maille_angle"] = "90.000000 90.000000 90.000000 deg"
                atoms_types[i]["maille_orient"] = "0"
                atoms_types[i]["maille_ref"] = "1"
        self.atoms_types = atoms_types

    def get_FFpar(self, atypes):
        """Return list of force field parameters."""
        ffparms = []
        # bonds
        bonds_list = []
        for ia, ja in it.combinations_with_replacement(atypes, 2):
            if "DU" in (ia, ja):
                continue
            elif "HT" in (ia, ja) and "HM" in (ia, ja):
                continue

            elif "HT" == ia and "HT" == ja:
                continue

            elif "HM" == ia and "HM" == ja:
                continue

            elif "CT" in (ia, ja) and "HM" in (ia, ja):
                continue

            elif "CM" in (ia, ja) and "HT" in (ia, ja):
                continue

            try:
                ffparms.append(BONDS_par_OPLS[(ia, ja)])
            except KeyError:
                ffparms.append(BONDS_par_OPLS[(ja, ia)])

            bonds_list.append((ia, ja))

        for ia, ja, ka in ANGLES_par_OPLS:
            if ia in atypes and ja in atypes and ka in atypes:
                ffparms.append(ANGLES_par_OPLS[(ia, ja, ka)])

        for ia, ja, ka, la in DIHEDRALS_par_OPLS:
            if ia in atypes and ja in atypes and ka in atypes and la in atypes:
                ffparms.append(DIHEDRALS_par_OPLS[(ia, ja, ka, la)])

        self.ffparms = ffparms

    @decoTime
    def save_file(self):
        """Save the information in a file."""
        print("Writing FAtomes file:", end=" - ")
        lines = ""
        now = datetime.now()
        now = now.strftime("%d/%m/%y %H:%M:%S")
        header = f"""*
*#####################################################################*
*                       Fichier issu de Stamp                         *
*                        (Iteration   000000)                         *
*                                                                     *
*        Date: {now}                                      *
*        Created using STAMPtools                                     *
*#####################################################################*
*
* ======================= 
* Bibliotheque des atomes 
* =======================
*
"""
        lines += header

        # Atoms types
        # ==============================
        atoms_types = self.atoms_types
        lines += "NbTypesAtomes %d\n" % len(atoms_types)

        for i in atoms_types:
            type_lines = f"* Description du type atomique    {i + 1}\n"
            for key in atoms_types[i]:
                if atoms_types[i][key] != "":
                    type_lines += "{} {}\n".format(key, atoms_types[i][key])
            
            lines += type_lines

        # Intermoleculars potentials
        # ==============================
        lines += """*
* ============================ 
* Potentiels intermoleculaires 
* ============================ 
"""
        pairsPot = {}
        for i in atoms_types:
            atom = atoms_types[i]["nomFF"]
            pairsPot[atom] = {}
            try:
                pairsPot[atom][atom] = POTENTIAL_par_OPLS[atom]
            except KeyError:
                pairsPot[atom][atom] = POTENTIAL_par_GAFF[atom]

        for atom1 in pairsPot:
            for atom2 in pairsPot[atom1]:
                pars = pairsPot[atom1][atom2]
                lines += f"Potentiel  {atom1}  {atom2}  {pars}\n"

        # Forcefield parameters
        # ==============================
        lines += """*
* =============== 
* Champ de forces 
* =============== 
ChampDeForces
"""
        self.get_FFpar(pairsPot.keys())
        ffparms = self.ffparms
        lines += "%d\n" % len(ffparms)
        for p in ffparms:
            lines += "%s\n" % p

        # Atomics positions
        # ==============================
        lines += """*
* =================== 
* Positions atomiques 
* =================== 
PositionDesAtomesCart ang
"""
        dfatoms = self.bulk.dfatoms.copy()
        dfatoms["x"] *= 10
        dfatoms["y"] *= 10
        dfatoms["z"] *= 10
        lines += "%d\n" % len(dfatoms)
        for i in dfatoms.index:
            line_xyz = (dfatoms.name[i], dfatoms.x[i], dfatoms.y[i], dfatoms.z[i])
            lines += '%5s%15.6f%15.6f%15.6f\n' % line_xyz

        # System connectivity
        # ==============================
        lines += """*
* ============ 
* Connectivite 
* ============ 
Zmatrice
"""
        connect = self.bulk.connect
        lines += "%d\n" % len(connect)
        for i in connect:
            if dfatoms.loc[i, "name"] == "DU":
                lines += "{}\n".format(i)
            else:
                lines += "{} {}\n".format(i, " ".join([str(j) for j in connect[i]]))

        # System reactivity
        # ==============================
        lines += """*
* ============ 
* Reactivite 
* ============ 
Reactivite
* INERTE/ACCEPTEUR/DONNEUR
"""
        atomsReact = dfatoms[dfatoms["name"] == "CT"].copy()
        atomsReact["reactivity"] = ["ACCEPTEUR", "DONNEUR"] * int(len(atomsReact) / 2)

        lines += "%d\n" % len(atomsReact)
        for i in atomsReact.index:
            lines += "{} {}\n".format(i, atomsReact.loc[i, "reactivity"])

        # End
        # ==============================
        with open(self.name, "w") as out:
            out.write(lines)


class TOPOL:
    """Object describing the system topology."""

    def __init__(self, fatomes):
        """Initialize the class by reading the indicated Fatomes file."""
        self.file = fatomes
        self._read_fatomes()

        # searching conectivity
        self._get_connectivity()

        # For remove DU
        self.key_remove = None

        # Add PC is actif
        self.adding_PC_proc = False

        # Add -OH terminal actif
        self.adding_OH_proc = False

    @decoTime
    def _read_fatomes(self):
        """Read Fatomes file."""
        print("Reading fatomes file", end=" - ")

        natypes = 0
        atomsM = {}
        xyz = []
        lnom = []
        lmass = []
        lcharges = []
        lcharge_mod = []
        connects = dict()

        # fatome_lines = []

        # ATOMS types
        atoms_types = {}
        infoType = {
            "nom": "",
            "nomXYZ": "",
            "nomFF": "",
            "type": "",
            "masse": "",
            "charge": "",
            "gele": "",
            "structure": "",
            "maille_long": "",
            "maille_angle": "",
            "maille_orient": "",
            "maille_ref": ""
        }
        ntype = 0
        natoms = 0

        # Intermolecular parameters
        Intermol_potentials = []
        comb_rules = ""

        # Force field parameters
        ffparms = []

        # Charges modified
        ChargesMOD = []

        # Contributions intramolecular to dispersion
        ContribDispRepIntra = []

        with open(self.file, "r") as FATM:
            nom = ""
            for line in FATM:
                # fatome_lines.append(line)
                if line.startswith("*"):
                    # ignore lines with the * symbol
                    continue

                if noms.match(line):
                    nom = ""
                    m = noms.match(line)
                    lnom.append(m.groupdict())
                    nom += m.groupdict()["nom"]
                    ntype += 1

                if ntype > 0 and ntype <= natypes:
                    break_line = line.split()
                    keywords = break_line[0]
                    if keywords in infoType:
                        atoms_types[ntype - 1][keywords] = " ".join(break_line[1:])

                if masses.match(line):
                    m = masses.match(line)
                    lmass.append(m.groupdict())

                if charges.match(line):
                    m = charges.match(line)
                    lcharges.append((nom, float(m.groupdict()["charge"])))
                    nom = ""

                if "NbTypesAtomes" in line:
                    line = line.split()
                    natypes += int(line[1])
                    atoms_types = {idx: dict(infoType) for idx in range(0, natypes)}
                    continue

                if "maille_long" in line:
                    line = line.split()
                    box = np.array(line[1:4]).astype(np.float64)
                    continue

                if atoms.match(line):
                    m = atoms.match(line)
                    xyz.append(m.groupdict())

                if "ChampDeForces" in line:
                    N = int(FATM.readline())
                    for _ in range(N):
                        lne = FATM.readline()
                        while lne.startswith("*"):
                            lne = FATM.readline()

                        ffparms.append(lne.replace("\n", ""))

                if "PositionDesAtomesCart" in line:
                    natoms += int(FATM.readline())

                if "Zmatrice" in line:
                    N = int(FATM.readline())
                    # print("N conectivity:", N)
                    for _ in range(N):
                        zline = FATM.readline()
                        zline = zline.split()
                        zline = [int(i) for i in zline]
                        connects[zline[0]] = zline[1:]
                if "VitesseDesAtomes" in line:
                    for _ in range(natoms):
                        FATM.readline()

                if "Regle_melange" in line:
                    comb_rules = line.replace("\n", "").split()[1]

                if "ModificationChargeDesAtomes" in line:
                    N = int(FATM.readline())
                    for _ in range(N):
                        chline = FATM.readline()
                        ChargesMOD.append(chline.replace("\n", ""))
                        chline = chline.split()
                        lcharge_mod.append({"idx": int(chline[0]), "charge": float(chline[1])})

                if "ContribDispRepIntra" in line:
                    N = int(FATM.readline())
                    for _ in range(N):
                        lne = FATM.readline()
                        while lne.startswith("*"):
                            lne = FATM.readline()

                        ContribDispRepIntra.append(lne.replace("\n", ""))

                if line.startswith("Potentiel"):
                    Intermol_potentials.append(line.replace("Potentiel  ", "").replace("\n", ""))

        lcharges = {key: charge for (key, charge) in lcharges}

        def from_lcharges(x):
            """Return partial charge for the atomic index."""
            try:
                return lcharges[x]
            except KeyError:
                return 0.0

        lcharge_mod = pd.DataFrame(lcharge_mod)
        atomsM = dict()
        if len(lmass) == len(lnom):
            for i in range(len(lnom)):
                atomsM[lnom[i]["nom"]] = np.float64(lmass[i]["mass"])
                if lnom[i]["nom"][0].upper() == "H":
                    atomsM[lnom[i]["nom"]] = np.float64(1.008e-03)
        else:
            print("\nERROR, no; and masses dont similar")
            exit()

        tabXYZ = pd.DataFrame(xyz)
        tabXYZ["atsb"] = tabXYZ["atypes"].apply(change_atsb)
        tabXYZ = tabXYZ.astype({
            "x": np.float64,
            "y": np.float64,
            "z": np.float64
        })

        tabXYZ["mass"] = 0.0
        for i in tabXYZ.index:
            try:
                tabXYZ["mass"] = atomsM[tabXYZ.loc[i, "atypes"]]
            except KeyError:
                print("ERROR")
                print("Atom type no reconize:", tabXYZ.loc[i, "atypes"])
                exit()

        tabXYZ["charge"] = tabXYZ["atypes"].apply(from_lcharges)
        if len(lcharge_mod) <= len(tabXYZ):
            for i in lcharge_mod.index:
                tabXYZ.loc[lcharge_mod.loc[i, "idx"], "charge"] = lcharge_mod.loc[i, "charge"]

        self.dfatoms = tabXYZ.copy()
        self.box = box
        self.connects = connects
        self.atoms_types = atoms_types
        pairsPot = {}
        for pot in Intermol_potentials:
            pot = pot.split()
            pars = "  ".join(pot[2:])
            try:
                # Using index
                iat = int(pot[0])
                jat = int(pot[1])

                itype = atoms_types[iat]["nomFF"]
                jtype = atoms_types[jat]["nomFF"]

                if (iat, itype) not in pairsPot:
                    pairsPot[(iat, itype)] = {}

                if (jat, jtype) not in pairsPot[(iat, itype)]:
                    pairsPot[(iat, itype)][(jat, jtype)] = pars

            except ValueError:
                # Using atom type
                itype = pot[0]
                jtype = pot[1]

                if itype not in pairsPot:
                    pairsPot[itype] = {}

                if jtype not in pairsPot[itype]:
                    pairsPot[itype][jtype] = pars

                pairsPot[itype][jtype] = pars

        self.Intermol_potentials = Intermol_potentials
        self.pairsPot = pairsPot
        self.comb_rules = comb_rules
        self.ffparms = ffparms
        self.ChargesMOD = ChargesMOD
        self.ContribDispRepIntra = ContribDispRepIntra

    @decoTime
    def _get_connectivity(self):
        """Brings system connectivity."""
        print("Searching connectivity", end=" - ")
        conn = structure.connectivity()
        conn.define_atoms(self.dfatoms)
        conn.read_dict(self.connects)
        self.connectivity = conn

    @property
    def atoms_per_mol(self):
        """List atoms per molecule."""
        return self.connectivity.atomsMOL

    def export_xyz(self, file="AtomesCart"):
        """Extract and save the atomic coordinates from the FAtomes."""
        coord = self.dfatoms.copy()
        coord["atsb"] = coord["atsb"].apply(change_atsb)
        structure.save_xyz(coord, file)
        print(f"Saved file: {file}")

    @decoTime
    def noPBC(self, xyz=None, name="system_noPBC"):
        """Return a system without periodic boundary conditions."""
        print("No PBC system", end=" - ")
        atoms_per_mol = self.atoms_per_mol
        connectivity = self.connectivity

        ###
        newconnectivity = structure.connectivity()
        conn_list = []
        for resid in atoms_per_mol:
            atoms = atoms_per_mol[resid]["index"]

            # Extracts the connectivity of the residue.
            mol_conn = connectivity.sub_connect(atoms)
            mol_conn.simple_at_symbols()

            # Reset index and symbols, and add mass
            mol_conn = mol_conn.reset_nodes()

            # If the residue has only one atom, it continues to the following.
            if len(atoms) == 1:
                conn_list.append(mol_conn)
                continue

            # whole molecule
            mol_conn.noPBC(self.box, center=np.zeros(3))

            # ncoords.append(mol_xyz)
            conn_list.append(mol_conn)

        # Create a new connectivity object with the new information.
        for mol in conn_list:
            newconnectivity.add_residue(mol)

        ncoords = newconnectivity.get_df()
        structure.save_xyz(ncoords, name=name)
        print(f"File save: {name}.xyz", end=" - ")

    @decoTime
    def complete_with_H(self):
        """Verify and complete molecules missing terminal hydrogens."""
        print("Adding H atom", end=" - ")
        atoms_per_mol = self.atoms_per_mol
        connectivity = self.connectivity

        ncoords = []
        newconnectivity = structure.connectivity()
        conn_list = []
        for resid in atoms_per_mol:
            atoms = atoms_per_mol[resid]["index"]
            # print("resid:", resid)

            # Extracts the connectivity of the residue.
            mol_conn = connectivity.sub_connect(atoms)
            mol_conn.simple_at_symbols()

            # Reset index and symbols, and add mass
            mol_conn = mol_conn.reset_nodes()

            # If the residue has only one atom, it continues to the following.
            if len(atoms) == 1:
                # mol_xyz = mol_conn.get_df()
                # ncoords.append(mol_xyz)
                conn_list.append(mol_conn)
                continue

            # whole molecule
            mol_conn.noPBC(self.box, center=np.zeros(3))

            # Search and add hydrogen to vacant atoms
            mol_conn.add_hydrogen(type_add="terminal-C", mass=1.008e-03, atypes="HT", charge=0.060)

            # Adding PBC
            mol_conn.addPBC(self.box, center=np.zeros(3))
            # mol_xyz = mol_conn.get_df()

            # ncoords.append(mol_xyz)
            conn_list.append(mol_conn)

        # Create a new connectivity object with the new information.
        for mol in conn_list:
            newconnectivity.add_residue(mol)

        self.connectivity = newconnectivity
        ncoords = self.connectivity.get_df()

        self.dfatoms = ncoords.copy()

    @decoTime
    def activate_PC(self, PCtopol):
        """Activate PC parameters in the FAtomes."""
        print("Activate PC", end=" - ")
        print("Topology selected:", PCtopol, end=" - ")
        dfatoms, connect = structure.load_mol2(PCtopol, connect=True)
        natomsPC = len(dfatoms)
        # Original dfatoms
        origDfatoms = self.dfatoms.copy()

        PCatoms = origDfatoms[origDfatoms["atypes"] == "DU"].copy()

        # How many PCs are in the system?
        nmolPC = int(len(PCatoms) / natomsPC)
        print("N PC in the system:", nmolPC, end=" - ")
        # Assigning atom types
        atypes = list(dfatoms["atsb"].values)
        PCatoms["atypes"] = atypes * nmolPC
        # Only for PC azo Orlando
        PCatoms["atsb"] = PCatoms["atypes"].apply(lambda x: x[0].upper())
        PCatoms["charge"] = list(dfatoms["charge"].values) * nmolPC
        PCatoms["mass"] = PCatoms["atsb"].apply(lambda at: structure.Elements[at]["mass"])

        toAddtypes = PCatoms["atypes"].unique()

        # Remove DU from atoms types
        atoms_types = self.atoms_types.copy()
        key_remove = None
        for i in atoms_types:
            if atoms_types[i]["nom"] == "DU":
                key_remove = i

        if key_remove is not None:
            del atoms_types[key_remove]
            self.key_remove = key_remove
        news_atoms_types = {}
        for i, atp in enumerate(toAddtypes):
            news_atoms_types[i] = {}
            news_atoms_types[i]["nom"] = atp
            try:
                news_atoms_types[i].update(ATOMSTYPES_par[atp])
            except KeyError:
                print("There is a type of atom that does not exist in the database, you must add it.")
                print("atom type:", atp)
                exit()

            if i == 0:
                news_atoms_types[i]["structure"] = "FICHIER"
                news_atoms_types[i]["maille_long"] = "{:.6f} {:.6f} {:.6f} ang".format(*tuple(self.box))
                news_atoms_types[i]["maille_angle"] = "90.000000 90.000000 90.000000 deg"
                news_atoms_types[i]["maille_orient"] = "0"
                news_atoms_types[i]["maille_ref"] = "1"

        for i, key in enumerate(atoms_types):
            news_atoms_types[i+len(toAddtypes)] = {}
            news_atoms_types[i+len(toAddtypes)] = atoms_types[key]

        # PC connect
        Npc_connects_map = {}
        for n_pc in range(nmolPC):
            for at in connect:
                Npc_connects_map[at + n_pc*len(connect)] = [c_at + n_pc*len(connect) for c_at in connect[at]]

        self.atoms_types = news_atoms_types
        self.dfatoms[self.dfatoms["atypes"] == "DU"] = PCatoms
        self.adding_PC_proc = True
        self.NPC_connects_map = Npc_connects_map
        self.ChargesMOD = list(self.dfatoms.loc[:len(Npc_connects_map)-1, "charge"])
        self.indexPC = slice(0, len(Npc_connects_map)-1)

    @decoTime
    def complete_with_OH(self):
        """Verify and complete molecules missing terminal hydrogens."""
        print("Adding OH atom", end=" - ")
        atoms_per_mol = self.atoms_per_mol
        connectivity = self.connectivity

        ncoords = []
        newconnectivity = structure.connectivity()
        conn_list = []

        for resid in atoms_per_mol:
            atoms = atoms_per_mol[resid]["index"]

            # Extracts the connectivity of the residue.
            mol_conn = connectivity.sub_connect(atoms)
            mol_conn.simple_at_symbols()

            # Reset index and symbols, and add mass
            mol_conn = mol_conn.reset_nodes()

            # If the residue has only one atom, it continues to the following.
            if len(atoms) == 1:
                # mol_xyz = mol_conn.get_df()
                # ncoords.append(mol_xyz)
                conn_list.append(mol_conn)
                continue

            # whole molecule
            mol_conn.noPBC(self.box, center=np.zeros(3))

            # Search and add oxygens to vacant carbons
            mol_conn.add_OH(type_add="terminal-C", mass=1.600e-02, atypes="OHopls", charge=-0.683)
            # Search and add hydrogen to vacant atoms
            mol_conn.add_hydrogen(type_add="OH", mass=1.008e-03, atypes="HOopls", charge=0.418)
            # mol_xyz = mol_conn.get_df()
            # structure.save_xyz(mol_xyz, "testOH")
            # Adding PBC
            mol_conn.addPBC(self.box, center=np.zeros(3))

            conn_list.append(mol_conn)
        
        # Create a new connectivity object with the new information.
        for mol in conn_list:
            newconnectivity.add_residue(mol)

        self.connectivity = newconnectivity
        ncoords = self.connectivity.get_df()
        if self.adding_PC_proc:
            ncoords.loc[self.indexPC, :] = self.dfatoms.loc[self.indexPC, :]

        self.dfatoms = ncoords.copy()
        self.adding_OH_proc = True

    @decoTime
    def center_to(self, ref, ref_type="atoms_ndx"):
        """Centers the system using a reference."""
        print("Centering the system to a reference", end=" - ")
        dfatoms = self.dfatoms.copy()
        if ref_type == "atoms_ndx":
            print("Reference is atoms index", end=" - ")
            dfref = dfatoms.loc[ref, :]
            ref_xyz = dfref.loc[:, ["x", "y", "z"]].values
            center = center_of_mass_polar(
                ref_xyz,
                self.box[0:3],
                np.ones(ref_xyz.shape[0]))

            dfatoms.loc[:, ["x", "y", "z"]] = translate_to(
                dfatoms.loc[:, ["x", "y", "z"]].values,
                center,
                self.box
            )

            dfatoms["x"] = dfatoms["x"] - center[0] + 0.5 * self.box[0]
            dfatoms["y"] = dfatoms["y"] - center[1] + 0.5 * self.box[1]
            dfatoms["z"] = dfatoms["z"] - center[2] + 0.5 * self.box[2]

            # update coordinates
            self.connectivity.update_coordinates(dfatoms)
            self.connectivity.addPBC(self.box, center=np.zeros(3))
            self.dfatoms = self.connectivity.get_df()

    def _verify_atomTypes(self):
        """Check the atom types with the modified atoms when adding new atoms."""
        new_types = []
        atoms_types = self.atoms_types.copy()
        NatomsTypes = len(atoms_types)
        modified_atoms = self.connectivity.modified_atoms
        if len(modified_atoms) > 0:
            for at in modified_atoms:
                if modified_atoms[at] not in new_types:
                    new_types.append(modified_atoms[at])

            for i, atyp in enumerate(new_types):
                atoms_types[i + NatomsTypes] = atyp

            self.atoms_types = atoms_types
            self.new_types = new_types

    def get_atomstypesPot(self):
        """Check that each type of atom contains an intermolecular van der waals potential."""
        atoms_types = self.atoms_types

        if self.adding_PC_proc:
            pairsPot = {}
        elif self.adding_OH_proc:
            pairsPot = {}
        else:
            pairsPot = self.pairsPot

        # key_remove = self.key_remove
        # if key_remove is not None:
        #     try:
        #         del pairsPot[key_remove]
        #     except KeyError:
        #         del pairsPot["DU"]
        #     finally:
        #         pass

        for atom in atoms_types:
            patom = (atom, atoms_types[atom]["nomFF"])
            if patom not in pairsPot:
                iatom = patom[0]
                itype = patom[1]
                if itype in POTENTIAL_par_GAFF:
                    pairsPot[patom] = {}
                    pairsPot[patom][patom] = POTENTIAL_par_GAFF[itype]

                elif itype in POTENTIAL_par_OPLS:
                    pairsPot[patom] = {}
                    pairsPot[patom][patom] = POTENTIAL_par_OPLS[itype]

                else:
                    print("Parameter not found:", iatom, itype)
                    exit()

        # Find the cross terms by applying the geometric combination rule.
        for (ati, typi), (atj, typj) in it.combinations(pairsPot, 2):
            if typi in NAMES_OPLS and typj in NAMES_OPLS:
                if (atj, typj) not in pairsPot[(ati, typi)]:
                    lj_pari = pairsPot[(ati, typi)][(ati, typi)]
                    lj_parj = pairsPot[(atj, typj)][(atj, typj)]

                    newPar = combRules_geom(
                        (float(lj_pari.split()[2]), float(lj_pari.split()[5])),
                        (float(lj_parj.split()[2]), float(lj_parj.split()[5]))
                    )

                    pairsPot[(ati, typi)][(atj, typj)] = "LJ  sigma  %.3f  ang  epsilon  %.3f  kcal/mol  rc  12.5  ang" % newPar

        self.pairsPot = pairsPot

    def check_newsFFpar(self):
        """List the force field parameters from the newly added parameters."""
        if self.adding_PC_proc:
            # angles
            for bond in BONDS_par_GAFF:
                self.ffparms.append(BONDS_par_GAFF[bond])
            
            # angles
            for ang in ANGLES_par_GAFF:
                self.ffparms.append(ANGLES_par_GAFF[ang])

            # torsion
            for dih in DIHEDRALS_par_GAFF:
                self.ffparms.append(DIHEDRALS_par_GAFF[dih])

            # impropers
            for imp in IMPROPERS_par_GAFF:
                self.ffparms.append(IMPROPERS_par_GAFF[imp])

        if self.adding_OH_proc:
            atypes = []
            pairsPot = self.pairsPot
            for iat, itype in pairsPot:
                atypes.append(itype)
            types_added = []
            for par in self.ffparms:
                if "bond_opls" in par:
                    types_added.append(tuple(par.split()[1:3]))
            new_types = [at["nomFF"] for at in self.new_types]
            for i, j in it.combinations(new_types, 2):
                if (i, j) in types_added or (j, i) in types_added:
                    pass
                elif i[0].upper() == "H" and j[0].upper() == "H":
                    pass
                elif i == "HT" and j == "OHopls":
                    pass
                elif j == "HT" and i == "OHopls":
                    pass
                elif i == "CT" and j == "HOopls":
                    pass
                elif j == "CT" and i == "HOopls":
                    pass
                else:
                    try:
                        self.ffparms.append(BONDS_par_OPLS[(i, j)])
                        types_added.append((i, j))
                    except KeyError:
                        self.ffparms.append(BONDS_par_OPLS[(j, i)])
                        types_added.append((j, i))

                    # angles
                    for ang in ANGLES_par_OPLS:
                        if i in ang or j in ang:
                            if ANGLES_par_OPLS[ang] not in self.ffparms:
                                self.ffparms.append(ANGLES_par_OPLS[ang])

                    # torsion
                    for dih in DIHEDRALS_par_OPLS:
                        if i in dih or j in dih:
                            if DIHEDRALS_par_OPLS[dih] not in self.ffparms:
                                self.ffparms.append(DIHEDRALS_par_OPLS[dih])

    @decoTime
    def save_fatomes(self, name=None):
        """Save the fatomes in a file."""
        print("Writing FAtomes file:", end=" - ")
        now = datetime.now()
        now = now.strftime("%d/%m/%y %H:%M:%S")
        header = f"""*
*#####################################################################*
*                       Fichier issu de Stamp                         *
*                        (Iteration   000000)                         *
*                                                                     *
*        Date: {now}                                      *
*        Created using STAMPtools                                     *
*#####################################################################*
*
* ======================= 
* Bibliotheque des atomes 
* =======================
*
"""

        out = ""
        if name:
            out += name.replace(".in", "") + ".in"
        else:
            out += self.file.replace(".in", "") + "_modif.in"

        lines = ""
        # header
        lines += header

        self._verify_atomTypes()
        atoms_types = self.atoms_types
        # Atoms types
        # ==============================
        lines += "NbTypesAtomes %d\n" % len(atoms_types)

        for i in atoms_types:
            type_lines = f"* Description du type atomique    {i}\n"
            for key in atoms_types[i]:
                if atoms_types[i][key] != "":
                    type_lines += "{} {}\n".format(key, atoms_types[i][key])
            
            lines += type_lines

        # Intermoleculars potentials
        # ==============================
        lines += """*
* ============================ 
* Potentiels intermoleculaires 
* ============================ 
"""
        self.get_atomstypesPot()
        pairsAdded = []
        pairsPot = self.pairsPot
        print(pairsPot)
        ######Cleaning new version
        ###newPairsPot = {}
        ###for patom1 in pairsPot:
        ###    iat = patom1[1]
        ###    newPairsPot[iat] = {}
        ###    for patom2 in pairsPot[patom1]:
        ###        pars = pairsPot[patom1][patom2]
        ###        jat = patom2[1]
        ###        newPairsPot[iat][jat] = pars
        ######
        ###for atp1 in newPairsPot:
        ###    for atp2 in newPairsPot[atp1]:
        ###        pars = newPairsPot[atp1][atp2]
        ###        lines += f"Potentiel  {atp1}  {atp2}  {pars}\n"

        ###
        for patom1 in pairsPot:
            for patom2 in pairsPot[patom1]:
                pars = pairsPot[patom1][patom2]
                iat = patom1[0]
                jat = patom2[0]
                if (iat, jat) not in pairsAdded or (jat, iat) not in pairsAdded:
                    lines += f"Potentiel  {iat}  {jat}  {pars}\n"
                    pairsAdded.append((iat, jat))

        if self.comb_rules != "":
            lines += "Regle_melange %s\n" % self.comb_rules
        elif self.adding_PC_proc:
            lines += "Regle_melange Lorentz-Berthelot\n"

        # Forcefield parameters
        # ==============================
        lines += """*
* =============== 
* Champ de forces 
* =============== 
ChampDeForces
"""     
        self.check_newsFFpar()
        ffparms = self.ffparms
        lines += "%d\n" % len(ffparms)
        for p in ffparms:
            lines += "%s\n" % p

        # Atomics positions
        # ==============================
        lines += """*
* =================== 
* Positions atomiques 
* =================== 
PositionDesAtomesCart ang
"""
        dfatoms = self.dfatoms.copy()
        lines += "%d\n" % len(dfatoms)
        for i in dfatoms.index:
            line_xyz = (dfatoms.atypes[i], dfatoms.x[i], dfatoms.y[i], dfatoms.z[i])
            lines += '%5s%15.6f%15.6f%15.6f\n' % line_xyz

        # System connectivity
        # ==============================
        lines += """*
* ============ 
* Connectivite 
* ============ 
Zmatrice
"""
        atoms_map = self.connectivity.atoms_map
        if self.adding_PC_proc:
            atoms_map.update(self.NPC_connects_map)
        lines += "%d\n" % len(atoms_map)
        for i in atoms_map:
            lines += "{} {}\n".format(i, " ".join([str(j) for j in atoms_map[i]]))

        # Atomic charges
        # ==============================
        lines_ch_title = """*
* =================
* Charges atomiques 
* =================
ModificationChargeDesAtomes e-
"""     
        if self.adding_OH_proc and not self.adding_PC_proc:
            lines_ch = ""
        else:
            ChargesMOD = self.ChargesMOD
            lines_ch = "%d\n" % len(ChargesMOD)  # (len(ChargesMOD) + len(self.connectivity.modified_atoms))
            for i, le in enumerate(ChargesMOD):
                lines_ch += "%d   %s\n" % (i, le)

            # for at in self.connectivity.modified_atoms:
            #     lines += "%d%10.3f\n" % (at, self.connectivity.nodes[at]["charge"])

        if len(lines_ch) > 0:
            lines += lines_ch_title
            lines += lines_ch

        # Contribution to dispersion intramolecular
        # ==============================
        lines_contr_title = """*
* =========================================================
* Contribution de dispersion repulsion en intramoleculaire 
* =========================================================
ContribDispRepIntra
"""     
        if self.adding_PC_proc:
            lines_contr = ContribDispRepIntra_PC
        elif self.adding_OH_proc:
            lines_contr = ""
        else:
            ContribDispRepIntra = self.ContribDispRepIntra
            lines_contr = "%d\n" % len(ContribDispRepIntra)
            for le in ContribDispRepIntra:
                lines_contr += "%s\n" % le

        if len(lines_contr) > 0:
            lines += lines_contr_title
            lines += lines_contr

        # Save lines
        with open(out, "w") as FATM:
            FATM.write(lines)

        print(f"file name: {out}", end=" - ")
