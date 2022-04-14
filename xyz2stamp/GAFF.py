
from xyz2stamp import structure as material
import numpy as np
import pandas as pd
import re


class ForceFieldError(Exception):
    pass


def assinging_at_type(atom):
    m_error = """\033[1;31m

    The force field for {} - {} is not defined
    Bonded with {}
    \033[m""".format(atom.n, atom.sb, atom.atoms_connect)

    if atom.sb == "C":
        """ CARBON """
        if atom.atoms_connect.count("H") == 4:
            return "C3"

        else:
            raise ForceFieldError(m_error)

    elif atom.sb == "H":
        """ HYDROGEN """
        if atom.atoms_connect == "Csp3":
            return "HC"

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

    dftypes["sigma"] = 0
    dftypes["epsilon"] = 0

    dfbonds["b0"] = 0
    dfbonds["kb"] = 0
    dfbonds["types"] = dfbonds["types"].apply(lambda x: tuple(sorted(x)))

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
        (?P<b0>[+-]?\d+\.\d+)\s+              # b0
        (?P<ln_kb>[+-]?\d+\.\d+)\s+           # ln force ctte
        """, re.X)
    bondspar = {}

    # ANGLES parameters
    angles = re.compile(
        r"""
        ^(?P<a1>\w+)\s-\s(?P<a2>\w+)\s-\s(?P<a3>\w+)\s+
        (?P<th0>[+-]?\d+\.\d+)\s+
        """, re.X)
    anglespar = {}

    with open(ffdata, "r") as DBASE:
        for line in DBASE:
            if vdw.match(line):
                # VDW
                m = vdw.match(line)
                m = m.groupdict()
                # print(m)

                if ff == "gaff":
                    vdwpar[m["type"]] = [
                            np.float64(m["sigma"]), np.float64(m["epsilon"])
                        ]

            elif bonds.match(line):
                # BONDS
                m = bonds.match(line)
                m = m.groupdict()
                bond = (m["b1"], m["b2"])

                if ff == "gaff":
                    bondspar[bond] = [
                            np.float64(m["b0"]), np.exp(np.float64(m["ln_kb"]))
                        ]

            elif angles.match(line):
                # ANGLES
                m = angles.match(line)
                m = m.groupdict()
                angle = (m["a1"], m["a2"], m["a3"])

                if ff == "gaff":
                    anglespar[angle] = [
                            np.float64(m["th0"])
                        ]

                    anglespar[angle[::-1]] = [
                            np.float64(m["th0"])
                        ]

    # VDW
    dftypes["sigma"] = dftypes["type"].apply(lambda x: vdwpar[x][0])
    dftypes["epsilon"] = dftypes["type"].apply(lambda x: vdwpar[x][1])
    # BONDS
    dfbonds["b0"] = dfbonds["types"].apply(lambda x: bondspar[x][0])
    dfbonds["kb"] = dfbonds["types"].apply(lambda x: bondspar[x][1])
    for i in dfbonds.index:
        ai = dfbonds.loc[i, "list"][0]
        aj = dfbonds.loc[i, "list"][1]
        MOL.connect[ai][aj]["dist"] = dfbonds.loc[i, "b0"]
        MOL.connect[aj][ai]["dist"] = dfbonds.loc[i, "b0"]

    # ANGLES
    dfangles["th0"] = dfangles["types"].apply(lambda x: anglespar[x][0])
    gaff_Kth(MOL)


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
    coord["charge"] = 0

    # Central Iterator
    for i in coord.index:
        if coord.atsb[i] == "C":
            """ CARBON """
            C = material.ATOM("C", i)
            coord.loc[C.n, "type"] = assinging_at_type(C)

        elif coord.atsb[i] == "H":
            """ HYDROGEN """
            H = material.ATOM("H", i)
            coord.loc[H.n, "type"] = assinging_at_type(H)

        else:
            raise ForceFieldError("""\033[1;31m

                It was not possible to assign atom type to %s
                \033[m""" % coord.atsb[i])

    MOL.dftypes = coord

    FFcoef(MOL, ffdata, ff="gaff")


def get_bonds(MOL, ffdata, databonds):
    """
    Get parameters for bonds.

    """
    i = 1
    idbond = list()
    bond_list = list()
    bond_type = list()

    coord = MOL.dftypes

    for iat, jat in MOL.bonds_list:
        if iat < jat:
            if (coord.loc[iat, 'ndx'], coord.loc[jat, 'ndx']) in databonds or (coord.loc[jat, 'ndx'], coord.loc[iat, 'ndx']) in databonds:
                idbond.append(i)
                bond_list.append((iat, jat))
                tp1 = coord.loc[iat, 'type']
                tp2 = coord.loc[jat, 'type']
                # bond_type.append(
                #     (re.sub("CH\d", "CHn", tp1), re.sub("CH\d", "CHn", tp2))
                # )
                i += 1

    dfbonds = pd.DataFrame({
        'bdid': idbond,
        'list': bond_list,
        'type': bond_type}
    )
    dfbonds = dfbonds.set_index('bdid')
    print(dfbonds)
    # dfbonds["code"] = "No assigned"
    # dfbonds["b0"] = "No assigned"
    # dfbonds["kb"] = "No assigned"

    # return dfbonds
