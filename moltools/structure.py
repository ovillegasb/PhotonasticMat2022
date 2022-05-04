
import re
import pandas as pd
import numpy as np
import networkx as nx
from scipy.spatial.distance import cdist
from scipy.constants import e


class MoleculeDefintionError(Exception):
    pass


# types of files sopported
typesfiles = ["xyz"]

# Defining error messages
_errorMessages = {
    0: """
    \033[1;31mMolecule format is not recognized,
    files supported to date: {}
    \033[m""".format(", ".join(typesfiles)),
    1: """
    \033[1;31mNo structures have been loaded yet.
    \033[m""",
    2: """
    \033[1;31mATOM is not configured to work with the atom: {}
    \033[m""",
    3: """
    \033[1;31mHybridation for {} - {} not found"
    \033[m""",
}


Elements = {  # g/mol
    'H': {'mass': 1.0080, 'num': 1},
    'C': {'mass': 12.011, 'num': 6},
}


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

    @property
    def atoms_map(self):
        """
        Returns a dict with the indices as key and the values as the atom
        indices to which it is bonded.
        """
        return {at: list(self.neighbors(at)) for at in self.nodes}

    def reset_nodes(self):
        """Reset le count of node from 0 to sizes nodes -1"""
        mapping = {value: count for count, value in enumerate(self.nodes, start=0)}

        # return nx.relabel_nodes(self, mapping, copy=True)
        self = nx.relabel_nodes(self, mapping, copy=True)

    def nbonds(self, inode):
        """Return number of atoms in connected to iat."""

        return int(self.degree[inode] / 2)

    def at_connect(self, inode):
        """Return the atoms symbols conectec to i node"""

        return '-'.join([self.nodes[a]['atsb'] for a in list(self.neighbors(inode))])

    def get_df(self):
        """Return a coordinate dataframe from connectivity."""
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


class MOL:
    """
    Class used to represent a molecule.

    Attributes
    ----------
    dfatoms : DataFrame
        Table with all the information of the atoms of the molecule, symbol,
        mass and x and z coordinates.

    Methods
    -------

    """

    # Atoms
    dfatoms = pd.DataFrame()
    dftypes = pd.DataFrame()

    atoms_count = {"C": 0, "N": 0, "H": 0}

    dfbonds = pd.DataFrame()
    dfangles = pd.DataFrame()
    dfdih = pd.DataFrame()
    dfimp = pd.DataFrame()

    # Connectivity
    connect = connectivity()

    bonds_list = []
    angles_list = []
    dihedrals_list = []
    impropers_list = []

    def load_file(self, file):
        """Extract information from file."""

        if file.endswith("xyz"):
            MOL.dfatoms = _load_xyz(file)

        elif file.endswith("out") or file.endswith("log"):
            MOL.dfatoms = _load_gaus(file)

        else:
            raise MoleculeDefintionError(_errorMessages[0])

        # RES name from file name
        self.res = file.split('/')[-1].split('.')[0].upper()

        MOL.atoms_count.update(dict(MOL.dfatoms["atsb"].value_counts()))

    def search_connectivity(self):
        """Search connectivity from a dfatoms and stores it in an object
        variable connect."""
        if len(MOL.dfatoms) > 0:
            MOL.connect.get_connectivity(MOL.dfatoms)
        else:
            raise MoleculeDefintionError(_errorMessages[1])

    @property
    def get_info(self):
        """Returns a header with molecule information."""
        info = ''
        info += '#----------------------\n'
        info += '# Molecular info: {}\n'.format(self.res)
        info += '#----------------------\n'
        info += '# MM {:.3f} g/mol\n'.format(self.mm)
        info += '# DBE {:d}\n'.format(self.DBE)
        info += '# Formule {}\n'.format(self.formule)
        info += '# Dipolar {:.3f} D\n'.format(self.mu)
        info += '#----------------------\n'
        return info

    @property
    def mu(self):
        """Compute the moment dipolar."""
        x = MOL.dfatoms.x.values
        y = MOL.dfatoms.y.values
        z = MOL.dfatoms.z.values
        q = MOL.dfatoms.charge.values
        mu_x = np.sum(q * x * (4.8 / 1.6e-29) * e) * 1e-10
        mu_y = np.sum(q * y * (4.8 / 1.6e-29) * e) * 1e-10
        mu_z = np.sum(q * z * (4.8 / 1.6e-29) * e) * 1e-10
        return np.sqrt(mu_x**2 + mu_y**2 + mu_z**2)

    @property
    def DBE(self):
        """Return the DBE (Double Bond Equivalent)."""
        atoms = MOL.atoms_count
        return int((2 * atoms["C"] + atoms["N"] + 2 - atoms["H"]) / 2)

    @property
    def mm(self):
        """Return the Molar Mass."""
        return MOL.dfatoms['mass'].sum()

    @property
    def formule(self):
        """Return the molecular formule."""
        formule = ""
        atoms = MOL.atoms_count
        for at in sorted(atoms):
            if atoms[at] > 0:
                if atoms[at] == 1:
                    formule += at

                else:
                    formule += at
                    formule += str(atoms[at])

        return formule

    @property
    def totq(self):
        """Compute the total charge."""
        return self.dftypes["charge"].sum()


class ATOM(MOL):
    """
    Subclass to represent an atom of the parent class MOL.
    """
    def __init__(self, at, n):
        """Initialize with the atomic symbol and the index in the molecule."""
        self.atsb = at
        self.n = n

    @property
    def atoms_connect(self):
        """Returns a string with the symbols of the connected atoms and their
        hybridization."""
        hyb_atoms_connect = ""

        if self.atsb == "H":
            # HYDROGEN
            at_i = list(MOL.connect.neighbors(self.n))[0]
            sb = MOL.connect.nodes[at_i]["atsb"]
            nbonds = len(list(MOL.connect.neighbors(at_i)))

            if sb == 'C' and nbonds == 3:
                """ H, aromatic hydrogen """
                hyb_atoms_connect += "Csp2"

            elif sb == 'C' and nbonds == 4:
                """ H, aliphatic hydrogen """
                hyb_atoms_connect += "Csp3"

        elif self.atsb == "C":
            # CARBON
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
        else:
            raise MoleculeDefintionError(_errorMessages[2].format(sb))

        return hyb_atoms_connect

    @property
    def hyb(self):
        atms_hyb = {
            "C": {4: "sp3", 3: "sp2"}
        }
        try:
            return atms_hyb[self.atsb][len(self.connect[self.n])]
        except KeyError:
            raise MoleculeDefintionError(_errorMessages[3].format(
                self.n, self.atsb))


def _remove_coord_negative(coord):
    for p in ["x", "y", "z"]:
        if coord[p].min() < 0.0:
            coord[p] = coord[p] - coord[p].min()
    return coord


def _load_xyz(file):
    """Read a file xyz."""
    coord = pd.read_csv(
        file,
        sep=r'\s+',
        skiprows=2,
        header=None,
        names=['atsb', 'x', 'y', 'z'],
        dtype={'x': np.float64, 'y': np.float64, 'z': np.float64}
    )

    coord = _remove_coord_negative(coord)

    coord["mass"] = coord["atsb"].apply(lambda at: Elements[at]["mass"])
    coord["num"] = coord["atsb"].apply(lambda at: Elements[at]["num"])
    # This file has no partial charges .
    coord["charge"] = 0.0
    return coord


def _load_gaus(file):
    """Obtains the optimized geometry and partial charges of the system from
    output file gaussian."""

    bond = re.compile(r"""
            ^\s!\sR\d+\s+
            R\((?P<at1>\d+),(?P<at2>\d+)\)\s+   # bond list
            """, re.X)

    coord = re.compile(r"""
            \s+
            (?P<atid>\d+)\s+           # Atom id.
            (?P<atnum>\d+)\s+          # Atomic number.
            \d+\s+
            (?P<x>[+-]?\d+\.\d+)\s+    # Orthogonal coordinates for X.
            (?P<y>[+-]?\d+\.\d+)\s+    # Orthogonal coordinates for Y.
            (?P<z>[+-]?\d+\.\d+)\s+    # Orthogonal coordinates for Z.
            """, re.X)

    charge = re.compile(r"""
            \s+
            (?P<atid>\d+)\s+               # Atom id.
            (?P<atsb>\w+)\s+               # Atomic number.
            (?P<charge>[+-]?\d+\.\d+)\s+   # Partial charges
            """, re.X)

    IGaus = list()
    data1 = list()
    data2 = list()
    data3 = list()

    with open(file, 'r') as INPUT:
        for line in INPUT:
            IGaus.append(line)
    IGaus = dict(enumerate(IGaus))

    for n in IGaus:
        """ Parametre Optimized """
        if 'Optimized Parameters' in IGaus[n]:
            for a in IGaus:
                if a > n:
                    """ Bonds list"""
                    if bond.match(IGaus[a]):
                        m = bond.match(IGaus[a])
                        a1 = int(m.groupdict()['at1'])
                        a2 = int(m.groupdict()['at2'])
                        data2.append((a1, a2))

                    """ Coordinates """
                    if 'Standard orientation:' in IGaus[a]:
                        for b in IGaus:
                            if b > a:
                                if coord.match(IGaus[b]):
                                    m = coord.match(IGaus[b])
                                    data1.append(m.groupdict())

                    if 'ESP charges:' in IGaus[a]:
                        for c in IGaus:
                            if c > a + 1:
                                if charge.match(IGaus[c]):
                                    m = charge.match(IGaus[c])
                                    data3.append(m.groupdict())
                                else:
                                    break

    coord = pd.DataFrame(data1)
    atcharges = pd.DataFrame(data3)
    coord['atsb'] = atcharges['atsb']
    coord['charge'] = atcharges['charge']

    """ Formating """
    coord = coord.astype({
        'atid': np.int,
        'num': np.int,
        'x': np.float,
        'y': np.float,
        'z': np.float,
        'charge': np.float})
    coord = coord.set_index('atid')

    coord["mass"] = coord["atsb"].apply(lambda at: Elements[at]["mass"])

    """ Listing bond, angle """
    # bond_list = set(data2)

    """ connectivity"""
    # connect = get_connectivity(dfatoms, bond_list)

    """ angles """
    # angle_list = get_angles_list(len(dfatoms), bond_list)

    """ dihedrals and pairs """
    # dihedral_list, pair_list = get_dihedral_pairs_list(len(dfatoms), bond_list)

    """ impropers """
    # improper_list = get_improper_list(len(dfatoms), bond_list)

    # return dfatoms, bond_list, angle_list, dihedral_list, improper_list, pair_list, connect

    return coord
