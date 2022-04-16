
from xyz2stamp import structure as material
import numpy as np
import pandas as pd
import re


class ForceFieldError(Exception):
    pass


def assinging_at_type(atom):
    m_error = """\033[1;31m

    The force field for {} - {} is not defined
    Bonded with {}.
    \033[m""".format(atom.n, atom.sb, atom.atoms_connect)

    if atom.sb == "C":
        """ CARBON """
        if atom.atoms_connect.count("H") == 4:
            return "c3"

        elif atom.hyb == "sp3" and atom.atoms_connect.count("H") == 3:
            return "c3"

        elif atom.hyb == "sp3" and atom.atoms_connect.count("H") == 2:
            return "c3"

        else:
            raise ForceFieldError(m_error)

    elif atom.sb == "H":
        """ HYDROGEN """
        if atom.atoms_connect == "Csp3":
            return "hc"

    else:
        raise ForceFieldError(m_error)


"""
Parameters to estimate the bond angle bending force constantes in GAFF.

"""
ParEstAngles = {
    "H": {"C": None, "Z": 0.784},
    "C": {"C": 1.339, "Z": 1.183}
}


def gaff_Kth(MOL):
    dfangles = MOL.dfangles
    connect = MOL.connect

    dfangles["kth"] = 0
    for i in dfangles.index:
        ai = dfangles.loc[i, "list"][0]
        aj = dfangles.loc[i, "list"][1]
        ak = dfangles.loc[i, "list"][2]

        r_ij = connect[ai][aj]["dist"]
        r_jk = connect[aj][ak]["dist"]
        # rad
        th0 = dfangles.loc[i, "th0"] * np.pi / 180

        D = (r_ij - r_jk)**2 / (r_ij + r_jk)**2
        Zi = ParEstAngles[connect.nodes[ai]["atsb"]]["Z"]
        Cj = ParEstAngles[connect.nodes[aj]["atsb"]]["C"]
        Zk = ParEstAngles[connect.nodes[ak]["atsb"]]["Z"]
        # Kcal / mol / rad2
        kth = 143.9 * Zi * Cj * Zk * np.exp(-2 * D) / (r_ij + r_jk) / th0**2
        # REVISAR unidates
        dfangles.loc[i, "kth"] = kth  # * 4.184  # * 90**2 / np.pi**2


def FFcoef(MOL, ffdata, ff):

    dftypes = MOL.dftypes
    dfbonds = MOL.dfbonds
    dfangles = MOL.dfangles
    dfdih = MOL.dfdih

    # print(MOL.get_interctions_list())
    # MOL.get_interctions_list()

    dftypes["sigma"] = 0
    dftypes["epsilon"] = 0

    dfbonds["b0"] = 0
    dfbonds["kb"] = 0

    # VDW parameters
    vdw = re.compile(
        r"""
        ^(?P<type>\w+)\s+
        (?P<sigma>[+-]?\d+\.\d+e?[+-]?\d?\d?)\s+      # nm
        (?P<epsilon>[+-]?\d+\.\d+e?[+-]?\d?\d?)\s+   # kJ/mol
        """, re.X)
    vdwpar = {}

    # BONDS parameters
    bonds = re.compile(
        r"""
        ^(?P<b1>\w+)\s-\s(?P<b2>\w+)\s+
        (?P<kb>[+-]?\d+\.\d+)\s+           # kb kcal/mol/ang2
        (?P<b0>[+-]?\d+\.\d+)\s+           # b0 ang
        """, re.X)
    bondspar = {}

    # ANGLES parameters
    angles = re.compile(
        r"""
        ^(?P<a1>\w+)\s-\s(?P<a2>\w+)\s-\s(?P<a3>\w+)\s+
        (?P<kth>[+-]?\d+\.\d+)\s+          # kth kcal/mol/rad2
        (?P<th0>[+-]?\d+\.\d+)\s+          # degree
        """, re.X)
    anglespar = {}

    # DIHEDRALS parameters
    dihedrals = re.compile(
        r"""
        ^(?P<d1>\w+)\s-\s(?P<d2>\w+)\s-\s(?P<d3>\w+)\s-\s(?P<d4>\w+)\s+
        (?P<divider>[+-]?\d+)\s+               # int
        (?P<Vn>[+-]?\d+\.\d+)\s+               # kcal/mol
        (?P<phi>[+-]?\d+\.\d+)\s+              # degree
        (?P<n>[+-]?\d\.\d*)\s+                # periodicity
        """, re.X)

    dihedralspar = {}

    with open(ffdata, "r") as DBASE:
        for line in DBASE:
            if vdw.match(line):
                # VDW
                m = vdw.match(line)
                m = m.groupdict()
                # print(m)
                # exit()

                if ff == "gaff":
                    vdwpar[m["type"]] = [
                            np.float64(m["sigma"]), np.float64(m["epsilon"])
                        ]

            elif bonds.match(line):
                # BONDS
                m = bonds.match(line)
                m = m.groupdict()
                bond = (m["b1"], m["b2"])
                # print(m)
                # print(bond)
                # exit()

                if ff == "gaff":
                    bondspar[bond] = [
                            np.float64(m["b0"]), np.float64(m["kb"])
                        ]

            elif angles.match(line):
                # ANGLES
                m = angles.match(line)
                m = m.groupdict()
                angle = (m["a1"], m["a2"], m["a3"])
                # print(m)
                # print(angle)
                # exit()

                if ff == "gaff":
                    anglespar[angle] = [
                            np.float64(m["th0"]), np.float64(m["kth"])
                        ]

                    anglespar[angle[::-1]] = [
                            np.float64(m["th0"]), np.float64(m["kth"])
                        ]

            elif dihedrals.match(line):
                # DIHEDRALS
                m = dihedrals.match(line)
                m = m.groupdict()
                dih = (m["d1"], m["d2"], m["d3"], m["d4"])
                # print(m)
                # print("HOLA" * 80)
                # print(dih)
                # print(dih)
                # exit()
                if ff == "gaff":
                    dihedralspar[dih] = [
                            np.int64(m["divider"]),
                            np.float64(m["Vn"]),
                            np.float64(m["phi"]),
                            int(np.float64(m["n"]))
                        ]

                    dihedralspar[dih[::-1]] = [
                            np.int64(m["divider"]),
                            np.float64(m["Vn"]),
                            np.float64(m["phi"]),
                            int(np.float64(m["n"]))
                        ]
    print(dihedrals)
    # exit()
    # VDW
    dftypes["sigma"] = dftypes["type"].apply(lambda x: vdwpar[x][0])
    dftypes["epsilon"] = dftypes["type"].apply(lambda x: vdwpar[x][1])
    print(dftypes)
    # exit()
    # BONDS
    try:
        dfbonds["b0"] = dfbonds["types"].apply(lambda x: bondspar[x][0])
        dfbonds["kb"] = dfbonds["types"].apply(lambda x: bondspar[x][1])
    except KeyError:
        # Exist bonds not found
        bonds = set(dfbonds.types.values)
        not_found = []
        for b in bonds:
            if b not in bondspar:
                not_found.append(b)
        raise ForceFieldError("""\033[1;31m
                It was not possible to assign a bond type:{}
                \033[m""".format(not_found))

    for i in dfbonds.index:
        ai = dfbonds.loc[i, "list"][0]
        aj = dfbonds.loc[i, "list"][1]
        MOL.connect[ai][aj]["dist"] = dfbonds.loc[i, "b0"]
        MOL.connect[aj][ai]["dist"] = dfbonds.loc[i, "b0"]
    print(dfbonds)
    # exit()
    # ANGLES
    try:
        dfangles["th0"] = dfangles["types"].apply(lambda x: anglespar[x][0])
        dfangles["kth"] = dfangles["types"].apply(lambda x: anglespar[x][1])
    except KeyError:
        # Exist angles not found
        angles = set(dfangles.types.values)
        not_found = []
        for a in angles:
            if a not in bondspar:
                not_found.append(a)
        raise ForceFieldError("""\033[1;31m
                It was not possible to assign a angle type:{}
                \033[m""".format(not_found))
    # gaff_Kth(MOL)
    print(dfangles)
    # DIHEDRALS
    try:
        dfdih["divider"] = dfdih["types"].apply(lambda x: dihedralspar[x][0])
        dfdih["Vn"] = dfdih["types"].apply(lambda x: dihedralspar[x][1])
        dfdih["phi"] = dfdih["types"].apply(lambda x: dihedralspar[x][2])
        dfdih["n"] = dfdih["types"].apply(lambda x: dihedralspar[x][3])
    except KeyError:
        # Exist angles not found
        dih = set(dfdih.types.values)
        not_found = []
        for d in dih:
            if d not in dihedralspar:
                not_found.append(d)
        raise ForceFieldError("""\033[1;31m
                It was not possible to assign a angle type:{}
                \033[m""".format(not_found))
    print(dfdih)


def Get_ATypes(MOL, ffdata):
    """
    Generates Gaff atoms types from a atomic coordinates dataframe

    Parameters:
    -----------

        table : DataFrame
            Coordinates atomics

    """
    coord = MOL.dfatoms
    coord["type"] = "No Found"

    # Central Iterator
    for i in coord.index:
        if coord.atsb[i] == "C":
            """ CARBON """
            C = material.ATOM("C", i)
            coord.loc[C.n, "type"] = assinging_at_type(C)

        elif coord.atsb[i] == "H":
            """ HYDROGEN """
            H = material.ATOM("H", i)
            # print("Here")
            # print(H.n)
            # print(H.sb)
            coord.loc[H.n, "type"] = assinging_at_type(H)

        else:
            raise ForceFieldError("""\033[1;31m

                It was not possible to assign atom type to %s
                \033[m""" % coord.atsb[i])

    MOL.dftypes = coord
    # print(coord)
    # exit()

    MOL.dfbonds["types"] = MOL.dfbonds["list"].apply(
        lambda x: (
            MOL.dftypes.loc[x[0], "type"],
            MOL.dftypes.loc[x[1], "type"])
        )
    MOL.dfbonds["types"] = MOL.dfbonds["types"].apply(
        lambda x: tuple(sorted(x)))

    MOL.dfangles["types"] = MOL.dfangles["list"].apply(
        lambda x: (
            MOL.dftypes.loc[x[0], "type"],
            MOL.dftypes.loc[x[1], "type"],
            MOL.dftypes.loc[x[2], "type"])
        )

    MOL.dfdih["types"] = MOL.dfdih["list"].apply(
        lambda x: (
            MOL.dftypes.loc[x[0], "type"],
            MOL.dftypes.loc[x[1], "type"],
            MOL.dftypes.loc[x[2], "type"],
            MOL.dftypes.loc[x[3], "type"])
        )

    FFcoef(MOL, ffdata, ff="gaff")
