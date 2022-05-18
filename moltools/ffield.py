
import numpy as np
import networkx as nx
import os
import re
from moltools.structure import ATOM
from pandas import DataFrame


class ForceFieldError(Exception):
    pass


# types of files sopported
forcefields = ["gromos", "gaff", "oplsaa"]

_errorMessages = {
    0: """
    \033[1;31mForce field is not found,
    Force field supported to date: {}
    \033[m""".format(", ".join(forcefields)),
    1: """
    \033[1;31mThe force field for {} - {} is not defined
    Bonded with {}.
    \033[m""",
    2: """
    \033[1;31mIt was not possible to assign atom type to {}
    \033[m""",
    3: """
    \033[1;31mIt was not possible to assign atom type to :\n{}
    \033[m"""
}


# database force field path
_ffpath = os.path.dirname(os.path.realpath(__file__))


def ffdata(ff):
    """search the path of force field"""
    return os.path.join(_ffpath, "forcefields", "%s.dat" % ff)


def database(ff):
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
        ^(?P<b1>\w+)\s*-\s*(?P<b2>\w+)\s+
        (?P<kb>[+-]?\d+\.\d+)\s+           # kb kcal/mol/ang2
        (?P<b0>[+-]?\d+\.\d+)\s+           # b0 ang
        """, re.X)
    bondspar = {}

    # ANGLES parameters
    angles = re.compile(
        r"""
        ^(?P<a1>\w+)\s*-\s*(?P<a2>\w+)\s*-\s*(?P<a3>\w+)\s+
        (?P<kth>[+-]?\d+\.\d+)\s+          # kth kcal/mol/rad2
        (?P<th0>[+-]?\d+\.\d+)\s+          # degree
        """, re.X)
    anglespar = {}

    # DIHEDRALS parameters
    dihedrals = re.compile(
        r"""
        ^(?P<d1>\w+)\s*-\s*(?P<d2>\w+)\s*-\s*(?P<d3>\w+)\s*-\s*(?P<d4>\w+)\s+
        (?P<divider>[+-]?\d+)\s+               # int
        (?P<Vn>[+-]?\d+\.\d+)\s+               # kcal/mol
        (?P<phi>[+-]?\d+\.\d+)\s+              # degree
        (?P<n>[+-]?\d\.\d*)\s+                 # periodicity
        """, re.X)
    dihedralspar = {}

    # IMPROPERS parameters
    impropers = re.compile(
        r"""
        ^(?P<d1>\w+)\s*-\s*(?P<d2>\w+)\s*-\s*(?P<d3>\w+)\s*-\s*(?P<d4>\w+)\s+
        (?P<Vn>[+-]?\d+\.\d+)\s+                # kcal/mol
        (?P<phi>[+-]?\d+\.\d*)\s+               # degree
        (?P<n>[+-]?\d\.\d*)\s+                  # periodicity
        """, re.X)
    improperspar = {}
    # print(ffdata(ff))
    with open(ffdata(ff), "r") as DBASE:
        for line in DBASE:
            if vdw.match(line):
                # VDW
                m = vdw.match(line)
                m = m.groupdict()
                # print(m)
                # exit()
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

            elif impropers.match(line):
                # IMPROPERS
                m = impropers.match(line)
                m = m.groupdict()
                imp = (m["d1"], m["d2"], m["d3"], m["d4"])
                improperspar[imp] = [
                            np.float64(m["Vn"]),
                            np.float64(m["phi"]),
                            int(np.float64(m["n"]))
                        ]

    return {
        "vdw": vdwpar,
        "bonds": bondspar,
        "angles": anglespar,
        "dihedrals": dihedralspar,
        "impropers": improperspar
        }


def get_interactions_list(MOL):
    """Lists of interactions are generated."""

    connect = MOL.connect
    all_ps = dict(nx.algorithms.all_pairs_shortest_path_length(connect))
    all_paths = []

    for s in all_ps.keys():
        for e in all_ps[s].keys():
            if all_ps[s][e] == 1:
                all_paths += list(nx.algorithms.all_simple_paths(connect, s, e, cutoff=1))
            elif all_ps[s][e] == 2:
                all_paths += list(nx.algorithms.all_simple_paths(connect, s, e, cutoff=2))
            elif all_ps[s][e] == 3:
                all_paths += list(nx.algorithms.all_simple_paths(connect, s, e, cutoff=3))
    # BONDS, list, types
    bonds_list = [tuple(p) for p in all_paths if len(set(p)) == 2]
    for iat, jat in bonds_list:
        if iat < jat:
            MOL.bonds_list.append((iat, jat))
    MOL.dfbonds["list"] = MOL.bonds_list

    MOL.dfbonds["types"] = MOL.dfbonds["list"].apply(
        lambda x: (
            MOL.dftypes.loc[x[0], "type"],
            MOL.dftypes.loc[x[1], "type"])
        )
    MOL.dfbonds["types"] = MOL.dfbonds["types"].apply(
        lambda x: tuple(sorted(x)))
    # ANGLES, list, types
    angles_list = [tuple(p) for p in all_paths if len(set(p)) == 3]
    for iat, jat, kat in angles_list:
        if iat < kat:
            MOL.angles_list.append((iat, jat, kat))
    MOL.dfangles["list"] = MOL.angles_list

    MOL.dfangles["types"] = MOL.dfangles["list"].apply(
        lambda x: (
            MOL.dftypes.loc[x[0], "type"],
            MOL.dftypes.loc[x[1], "type"],
            MOL.dftypes.loc[x[2], "type"])
        )
    # DIHEDRALS, list, types
    dihedrals_list = [tuple(p) for p in all_paths if len(set(p)) == 4]
    for iat, jat, kat, lat in dihedrals_list:
        if iat < lat:
            # Remove dihedrals -ca-ca-
            # if MOL.dftypes.loc[jat, "type"] != "ca" and MOL.dftypes.loc[kat, "type"] != "ca":
            #    MOL.dihedrals_list.append((iat, jat, kat, lat))
            MOL.dihedrals_list.append((iat, jat, kat, lat))
    MOL.dfdih["list"] = MOL.dihedrals_list

    MOL.dfdih["types"] = MOL.dfdih["list"].apply(
        lambda x: (
            MOL.dftypes.loc[x[0], "type"],
            MOL.dftypes.loc[x[1], "type"],
            MOL.dftypes.loc[x[2], "type"],
            MOL.dftypes.loc[x[3], "type"])
        )

    # OUT-OF-PLANE, IMPROPER, list
    for at in MOL.dftypes.index:
        if MOL.dftypes.loc[at, "type"] in ["ca"]:
            i, j, z = set(nx.all_neighbors(MOL.connect, at))
            # print(i, j, z)
            test = MOL.dftypes.loc[[i, j, z], "type"]  # .sort_values(axis=1)

            test = DataFrame({"idx": test.index, "type": test.values})

            # print(test.set_index('Column_to_make_index'))
            # print()
            # print(test)
            test.sort_values("type", ascending=True, inplace=True)
            test.set_index("idx", inplace=True)
            i, j, z = tuple(test.index)

            # i
            # print(type(test))
            # print("#"*50)
            MOL.impropers_list.append((i, j, at, z))
    # re.sub("CH\d", "CHn", tp1)
    # exit()

    MOL.dfimp["list"] = MOL.impropers_list

    MOL.dfimp["types"] = MOL.dfimp["list"].apply(
        lambda x: (
            MOL.dftypes.loc[x[0], "type"],  # re.sub("c\w", "X", MOL.dftypes.loc[x[0], "type"]),
            MOL.dftypes.loc[x[1], "type"],  # re.sub("c\w", "X", MOL.dftypes.loc[x[1], "type"]),
            MOL.dftypes.loc[x[2], "type"],
            MOL.dftypes.loc[x[3], "type"])
        )

    # print(MOL.dfimp)
    # print("Estamos aqui")
    # exit()

    # print("Natoms:", connect.number_of_nodes())
    # print("Bonds:\n", MOL.bonds_list)
    # print("Nbonds:", len(MOL.bonds_list))
    # print("Angles:\n", MOL.angles_list)
    # print("Nangles:", len(MOL.angles_list))
    # print("Dihedrals:\n", MOL.dihedrals_list)
    # print("NDihedrals:", len(MOL.dihedrals_list))


def assinging_at_type(atom):
    error = _errorMessages[1].format(atom.n, atom.atsb, atom.atoms_connect)
    if atom.atsb == "C":
        """ CARBON """
        if atom.atoms_connect.count("H") == 4:
            return "c3"

        elif atom.hyb == "sp3" and atom.atoms_connect.count("H") == 3:
            return "c3"

        elif atom.hyb == "sp3" and atom.atoms_connect.count("H") == 2:
            return "c3"

        elif atom.hyb == "sp2" and atom.atoms_connect.count("Csp2") >= 2:
            return "ca"

        else:
            raise ForceFieldError(error)

    elif atom.atsb == "H":
        """ HYDROGEN """
        if atom.atoms_connect == "Csp3":
            return "hc"

        elif atom.atoms_connect == "Csp2":
            return "ha"

    elif atom.atsb == "N":
        """ NYTROGEN """
        if atom.hyb == "sp2" and atom.atoms_connect.count("Csp2") < 2:
            # o ne?
            return "n2"

    else:
        raise ForceFieldError(error)


def get_atoms_types(MOL, ff):
    """Assigns the atom timpos for the indicated force field."""
    if ff == "gromos":
        print("\033[1;35mForce field: GROMOS\033[m")
        # GROMOS.Get_ATypes(MOL, ffdata(ff))

    elif ff == "gaff":
        print("\033[1;35mForce field: GAFF\033[m")
        # ffdata = os.path.join(ffpath, "forcefields", "%s.dat" % self.ff)
        # GAFF.Get_ATypes(MOL, ffdata(ff))
        # GAFF.get_bonds(MOL, ffdata)
    elif ff == "oplsaa":
        print("\033[1;35mForce field: OPLS AS\033[m")

    else:
        raise ForceFieldError(_errorMessages[0])

    coord = MOL.dfatoms
    coord["type"] = "No Found"

    # Central Iterator
    for i in coord.index:
        if coord.atsb[i] == "C":
            """ CARBON """
            C = ATOM("C", i)
            coord.loc[C.n, "type"] = assinging_at_type(C)

        elif coord.atsb[i] == "H":
            """ HYDROGEN """
            H = ATOM("H", i)
            coord.loc[H.n, "type"] = assinging_at_type(H)

        elif coord.atsb[i] == "N":
            """ NYTROGEN """
            N = ATOM("N", i)
            coord.loc[N.n, "type"] = assinging_at_type(N)

        else:
            raise ForceFieldError(_errorMessages[2].format(coord.atsb[i]))

    # print(coord)
    MOL.dftypes = coord


def get_ffparameters(MOL, ff):
    dftypes = MOL.dftypes
    dfbonds = MOL.dfbonds
    dfangles = MOL.dfangles
    dfdih = MOL.dfdih
    dfimp = MOL.dfimp

    # VDW
    dbase = database(ff)
    vdwpar = dbase["vdw"]
    try:
        dftypes["sigma"] = dftypes["type"].apply(lambda x: vdwpar[x][0])
        dftypes["epsilon"] = dftypes["type"].apply(lambda x: vdwpar[x][1])
    except KeyError:
        raise ForceFieldError(_errorMessages[3].format(dftypes))
    print(dftypes)

    # BONDS
    bondspar = dbase["bonds"]
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

    # ANGLES
    anglespar = dbase["angles"]
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
    print(dfangles)

    # DIHEDRALS
    dihedralspar = dbase["dihedrals"]
    print(dfdih)
    exit()
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
                It was not possible to assign a dihedral type:{}
                \033[m""".format(not_found))
    print(dfdih)

    # IMPROPERS
    improperspar = dbase["impropers"]
    try:
        dfimp["Vn"] = dfimp["types"].apply(lambda x: improperspar[x][0])
        dfimp["phi"] = dfimp["types"].apply(lambda x: improperspar[x][1])
        dfimp["n"] = dfimp["types"].apply(lambda x: improperspar[x][2])
    except KeyError:
        # Exist angles not found
        imp = set(dfimp.types.values)
        not_found = []
        for d in imp:
            if d not in improperspar:
                not_found.append(d)
        raise ForceFieldError("""\033[1;31m
                It was not possible to assign a improper dihedral type:{}
                \033[m""".format(not_found))
    print(dfimp)
