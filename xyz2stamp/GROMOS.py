
# from .structure import MOL

from xyz2stamp import structure as material
import numpy as np
import pandas as pd


class ForceFieldError(Exception):
    pass


def assinging_at_type(atom):
    m_error = '\nThe force field for {} - {} is not defined\nBonded with {}'.format(atom.n, atom.sb, atom.atoms_connect)
    if atom.sb == "C":

        """ CARBON """
        if atom.atoms_connect.count("H") == 4:
            return "CH4"

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

    if "CH" in row["type"]:
        row["mass"] = massUA[row["type"]]
    row["ndx"] = i

    return row


def Get_ATypes():
    """
    Generates Gromos 54a7 atoms types from a atomic coordinates dataframe

    Parameters:
    -----------

        table : DataFrame
            Coordinates atomics

    """
    coord = material.MOL.dfatoms
    coord["type"] = "No Found"
    out = list()

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
                continue
            out.append(get_row(H, i))
            coord.loc[H.n, "type"] = "Found"

        else:
            raise ForceFieldError("It was not possible to assign atom type to %s" % coord.atsb[i])

    out = pd.DataFrame(out, index=np.arange(1, len(out) + 1))
    material.MOL.dftypes = out
