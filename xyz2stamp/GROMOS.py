
from xyz2stamp import structure as material
import numpy as np
import pandas as pd
import re


class ForceFieldError(Exception):
    pass


def assinging_at_type(atom):
    m_error = '\nThe force field for {} - {} is not defined\nBonded with {}'.format(atom.n, atom.sb, atom.atoms_connect)
    if atom.sb == "C":

        """ CARBON """
        if atom.atoms_connect.count("H") == 4:
            return "CH4"

        elif atom.atoms_connect.count("H") == 3:
            return "CH3"

        else:
            raise ForceFieldError(m_error)

    elif atom.sb == "H":
        """ HYDROGEN """
        if atom.atoms_connect == "Csp3":
            return None

    else:
        raise ForceFieldError(m_error)


def get_row(at, i):
    # contruye la fila del atomo con su typo, al final se uniran para formar
    # un solo dataframe

    massUA = {
        "CH4": 16.04300,
        "CH3": 15.03500,
        "CH2": 14.02700,
        "CH1": 13.01900
    }

    row = at.dfatoms.loc[at.n, :]
    row = row.to_dict()
    row["type"] = assinging_at_type(at)
    row["charge"] = at.charge

    # if row["type"]:
    try:
        if "CH" in row["type"]:
            row["mass"] = massUA[row["type"]]
        row["ndx"] = i
    except TypeError:
        print(row)
        exit()

    return row


def FFcoef(dftypes, ffdata, ff="gromos"):

    dftypes["sigma"] = 0
    dftypes["epsilon"] = 0
    """
    def get_atype(x):
        atype = m["type"]
        if atype == x:
            yield m["c6"]

            yield m["c12"]
    """

    def sigma(c6, c12):
        # nm --> angstrom
        sig = (c12/c6)**(1/6)
        return sig / 10

    def epsilon(c6, c12):
        # 1 kcal = 4.184 kj
        # kJ / mol --> kcal / mol
        return (c6**2) / (4 * c12)

    # VDW parameters
    vdw = re.compile(
        r"""
        ^(?P<type>\w+)\s+
        (?P<c6>[+-]?\d+\.\d+)\s+                    # kJ/mol nm6
        (?P<c12>[+-]?\d+\.\d+e?[+-]?\d?\d?)\s+      # kJ/mol nm12
        """, re.X)

    with open(ffdata, "r") as DBASE:
        for line in DBASE:
            # VDW
            if vdw.match(line):
                m = vdw.match(line)
                m = m.groupdict()
                print(m)

                if ff == "gromos":
                    m = {
                        m["type"]: [
                            sigma(np.float64(m["c6"]), np.float64(m["c12"])),  # sigma (nm)
                            epsilon(np.float64(m["c6"]), np.float64(m["c12"]))  # epsilon (Kj/mol)
                        ]
                    }
                # print(m)

                dftypes["sigma"] = dftypes["type"].apply(lambda x: m[x][0])
                dftypes["epsilon"] = dftypes["type"].apply(lambda x: m[x][1])

                # print(dftypes)


def get_vdw_par(MOL):
    pass


def Get_ATypes(MOL, ffdata):
    """
    Generates Gromos 54a7 atoms types from a atomic coordinates dataframe

    Parameters:
    -----------

        table : DataFrame
            Coordinates atomics

    """
    coord = MOL.dfatoms
    connect = MOL.connect
    print(connect)
    print(coord)
    coord["type"] = "No Found"
    out = list()

    # print(coord)
    # exit()
    to_remove = []

    # Central Iterator
    for i in coord.index:

        if coord.atsb[i] == "C":
            """ CARBON """
            C = material.ATOM("C", i)
            out.append(get_row(C, i))
            coord.loc[C.n, "type"] = "Found"

        elif coord.atsb[i] == "H":
            """ HYDROGEN """
            H = material.ATOM("H", i)

            if H.atoms_connect == "Csp3":
                coord.loc[H.n, "type"] = "Found"
                # connect.remove_node(i)
                to_remove.append(i)
                continue
            out.append(get_row(H, i))
            coord.loc[H.n, "type"] = "Found"

        else:
            raise ForceFieldError("It was not possible to assign atom type to %s" % coord.atsb[i])

    # out = pd.DataFrame(out, index=np.arange(1, len(out) + 1))
    out = pd.DataFrame(out)
    MOL.dftypes = out
    # print(coord)
    if len(to_remove) > 0:
        for i in to_remove:
            connect.remove_node(i)

    MOL.connect = connect

    # print(MOL.connect)
    # print(out)
    # exit()

    # FFcoef(MOL.dftypes, ffdata, ff="gromos")
