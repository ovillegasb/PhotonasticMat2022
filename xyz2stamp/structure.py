
import networkx as nx
import itertools as it
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import os

from xyz2stamp import GROMOS, GAFF


class MoleculeDefintionError(Exception):
    pass


# database force field path
ffpath = os.path.dirname(os.path.realpath(__file__))


class MOL:
    """Class to represent a molecule object."""
    dfatoms = pd.DataFrame()
    connect = None
    dftypes = pd.DataFrame()

    dfbonds = pd.DataFrame()
    dfangles = pd.DataFrame()
    dfdih = pd.DataFrame()

    bonds_list = []
    angles_list = []
    dihedrals_list = []

    def __init__(self, dfatoms, connect):
        """
        Initialize calling the dfatoms and connect.

        """
        MOL.dfatoms = dfatoms
        MOL.connect = connect

        self._get_interctions_list()

    def _get_interctions_list(self):
        """Lists of interactions are generated."""

        connect = self.connect
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

        bonds_list = [tuple(p) for p in all_paths if len(set(p)) == 2]
        for iat, jat in bonds_list:
            if iat < jat:
                MOL.bonds_list.append((iat, jat))

        MOL.dfbonds["list"] = MOL.bonds_list

        angles_list = [tuple(p) for p in all_paths if len(set(p)) == 3]
        for iat, jat, kat in angles_list:
            if iat < kat:
                MOL.angles_list.append((iat, jat, kat))

        MOL.dfangles["list"] = MOL.angles_list

        dihedrals_list = [tuple(p) for p in all_paths if len(set(p)) == 4]
        print(dihedrals_list)
        print(len(dihedrals_list))
        for iat, jat, kat, lat in dihedrals_list:
            if iat < lat:
                MOL.dihedrals_list.append((iat, jat, kat, lat))

        MOL.dfdih["list"] = MOL.dihedrals_list

        # pairs_list = [(p[0], p[3]) for p in all_paths if len(set(p)) == 4]

        print("Natoms:", connect.number_of_nodes())
        print("Bonds:\n", MOL.bonds_list)
        print("Nbonds:", len(MOL.bonds_list))
        print("Angles:\n", MOL.angles_list)
        print("Nangles:", len(MOL.angles_list))
        print("Dihedrals:\n", MOL.dihedrals_list)
        # for i in MOL.dihedrals_list:
        #     # print(MOL.dftypes.loc[list(i), :])
        #     print(i)
        #     print(MOL.dfatoms.loc[i, :])
        # exit()
        print("NDihedrals:", len(MOL.dihedrals_list))
        # print(pairs_list)


class ATOM(MOL):
    """
    Clase que contruye a un atomo, lo cual pertenece a una clase superior, o
    molecula.

    """

    def __init__(self, sb, n):
        self.sb = sb  # simbolo atomico
        self.n = n  # indice del atomo

    @property
    def atoms_connect(self):
        hyb_atoms_connect = ""

        if self.sb == "H":
            at_i = list(MOL.connect.neighbors(self.n))[0]
            sb = MOL.connect.nodes[at_i]["atsb"]
            nbonds = len(list(MOL.connect.neighbors(at_i)))

            if sb == 'C' and nbonds == 3:
                """ H, aromatic hydrogen """
                hyb_atoms_connect += "Csp2"

            elif sb == 'C' and nbonds == 4:
                """ H, aliphatic hydrogen """
                hyb_atoms_connect += "Csp3"

        elif self.sb == "C":
            # atoms bonded to carbon
            for at_i in MOL.connect.neighbors(self.n):
                sb = MOL.connect.nodes[at_i]["atsb"]
                nbonds = len(list(MOL.connect.neighbors(at_i)))

                if sb == "C" and nbonds == 3:
                    hyb_atoms_connect += "Csp2"

                elif sb == "C" and nbonds == 4:
                    hyb_atoms_connect += "Csp3"

                elif sb == "H":
                    hyb_atoms_connect = "H" + hyb_atoms_connect

        return hyb_atoms_connect

    @property
    def hyb(self):
        atms_hyb = {
            "C": {4: "sp3", 3: "sp2"}
        }
        try:
            return atms_hyb[self.sb][len(self.connect[self.n])]
        except KeyError:
            raise MoleculeDefintionError("Hybridation for {} - {} not found".format(self.n, self.sb))

    @property
    def charge(self):
        try:
            charge = self.dfatoms.charge[self.n]

            if self.dfatoms.atsb[self.n] == "C":
                if self.hyb == "sp3":
                    for c in self.connect[self.n]:
                        if self.dfatoms.atsb[c] == "H":
                            charge += self.dfatoms.charge[c]
                elif self.hyb == "sp2":
                    pass

                else:
                    raise MoleculeDefintionError("Hybridation for {} - {} is not defined".format(self.n, self.sb))

        except AttributeError:
            charge = 0.0

        return np.round(charge, decimals=2)


class connectivity(nx.DiGraph):
    """Building a class connectivity from directed graphs."""
    def __init__(self):
        super().__init__()

    def get_connectivity(self, coord):
        """Build connectivity from coordinates using nodes like atoms."""

        # Add nodes using atoms andsymbols and coordinates
        for i in coord.index:
            self.add_node(
                i,
                xyz=coord.loc[i, ['x', 'y', 'z']].values,
                atsb=coord.loc[i, 'atsb']
            )

        # Add edges like bonds
        # get pairs atoms bonded
        pairs, m = self._neighboring_pairs(coord)
        for i, j in pairs:
            self.add_edge(i, j)
            self.add_edge(j, i)

            self.edges[i, j]['dist'] = m[i, j]
            self.edges[j, i]['dist'] = m[i, j]

    def nbonds(self, inode):
        """Return number of atoms in connected to iat."""

        return int(self.degree[inode] / 2)

    def _neighboring_pairs(self, coord):
        """Return neighboring pairs"""

        xyz = coord.loc[:, ['x', 'y', 'z']].values.astype(np.float64)
        # atoms = xyz = coord.loc[:, ['atsb']].values
        # compute distance
        m = cdist(xyz, xyz, 'euclidean')
        m = np.triu(m)

        indexs = np.where((m > 0.) & (m <= 1.7))
        # print(indexs)
        # print(indexs[0])
        # print(coord["atsb"][indexs[0]])
        # print(coord["atsb"][indexs[1]])
        # exit()

        pairs = map(lambda in0, in1: (in0, in1), indexs[0], indexs[1])

        return pairs, m

    def at_connect(self, inode):
        """Return the atoms symbols conectec to i node"""

        return '-'.join([self.nodes[a]['atsb'] for a in list(self.neighbors(inode))])

    def get_df(self):
        """Return the connectivity as a Pandas DataFrame."""
        indexs = list(self.nodes)
        rows = list()

        for i in self.nodes:
            rows.append({
                'atsb': self.nodes[i]['atsb'],
                'x': self.nodes[i]['xyz'][0],
                'y': self.nodes[i]['xyz'][1],
                'z': self.nodes[i]['xyz'][2]
            })

        df = pd.DataFrame(rows, index=indexs)

        return df


class FField:
    """
    Objeto que construye los parametros del campo de fuerza.

    """

    def __init__(self, ff):
        """
        Se inicializa cargando los parametros del campo de fuerza y asignando
        los tipos de atomos.

        """

        print(ffpath)

        ffdata = os.path.join(ffpath, "forcefields", "%s.dat" % ff)

        print("Campo de fuerza seleccionado:", ff)
        print(ffdata)

        # with open(ffdata, "r") as F:
        #     for line in F:
        #         print(line)

        self.ffdata = ffdata
        self.ff = ff

    def get_atoms_types(self, MOL):

        if self.ff == "gromos":
            print("\033[1;35mForce field: GROMOS\033[m")
            ffdata = os.path.join(ffpath, "forcefields", "%s.dat" % self.ff)
            GROMOS.Get_ATypes(MOL, ffdata)

        elif self.ff == "gaff":
            print("\033[1;35mForce field: GAFF\033[m")
            ffdata = os.path.join(ffpath, "forcefields", "%s.dat" % self.ff)
            GAFF.Get_ATypes(MOL, ffdata)
            # GAFF.get_bonds(MOL, ffdata)
