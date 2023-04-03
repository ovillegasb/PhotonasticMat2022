
"""
Module that stores different functions and tools for working with STAMP files.

Orlando Villegas - 2022

"""

import re
import numpy as np
import pandas as pd
from pymatgen.core import Molecule

""" Regular expression that extracts matrix XYZ """
atoms = re.compile(r"""
        ^\s+
        (?P<atsb>[A-Za-z]+\d?\d?)\s+      # Atom name.
        (?P<x>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for X.
        (?P<y>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Y.
        (?P<z>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Z.
        """, re.X)


def read_fatomes(file):
    """Read Fatomes file."""
    natypes = 0
    atomsM = {}
    xyz = []
    connects = dict()
    with open(file, "r") as FATM:
        for line in FATM:
            if "*" == line[0]:
                # ignore lines with the * symbol
                continue

            elif "NbTypesAtomes" in line:
                line = line.split()
                natypes += int(line[1])
                continue

            elif "nom" in line:
                line = line.split()
                nom = line[1]
                continue

            elif "masse" in line:
                line = line.split()
                mass = np.float64(line[1])
                atomsM[nom] = mass
                continue

            elif "maille_long" in line:
                line = line.split()
                box = np.array(line[1:4]).astype(np.float64)
                continue

            elif "PositionDesAtomesCart" in line:
                Natoms = int(FATM.readline())
                continue

            elif atoms.match(line):
                m = atoms.match(line)
                xyz.append(m.groupdict())

            elif "Zmatrice" in line:
                N = int(FATM.readline())
                print("N conectivity:", N)
                for _ in range(N):
                    zline = FATM.readline()
                    zline = zline.split()
                    zline = [int(i) for i in zline]
                    connects[zline[0]] = zline[1:]

    print("Number of atoms in XYZ matrix:", Natoms)

    print("ATOMS types and Mass [Kg/mol]")
    print(atomsM)

    print("Box dimensions [angs]:")
    print(box)

    tabXYZ = pd.DataFrame(xyz)

    tabXYZ = tabXYZ.astype({
        "x": np.float64,
        "y": np.float64,
        "z": np.float64
    })

    tabXYZ["mass"] = tabXYZ["atsb"].apply(lambda x: atomsM[x])

    return tabXYZ, box, connects


def save_gro(table, box, name="coord.gro"):
    """Save coordinate to file *.gro from dataframe with x, y, z."""
    nat = len(table)
    gro = name

    GRO = open(gro, "w", encoding="utf-8")
    GRO.write("GRO FILE\n")
    GRO.write("%5d\n" % nat)
    for i in table.index:
        GRO.write("{:>8}{:>7}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
            str(table.loc[i, "resid"]) + table.loc[i, "resname"],
            table.loc[i, "atsb"],
            i + 1,
            table.loc[i, "x"] * 0.1,
            table.loc[i, "y"] * 0.1,
            table.loc[i, "z"] * 0.1)
        )
    GRO.write("   {:.5f}   {:.5f}   {:.5f}\n".format(
        box[0] * 0.1,
        box[1] * 0.1,
        box[2] * 0.1))

    GRO.close()
    print("\nSaved gro file: \033[1;36m%s\033[m writed\n" % gro)


def read_donnees(file):
    """Read DONNEES file."""
    parameters = {}
    with open(file, "r") as D:
        for line in D:
            if not line.startswith("*"):
                line = line.replace("\n", "")
                line = line.split()[0:2]
                parameters[line[0]] = line[1]

    return parameters


def clean_data(df):
    """Clean data by removing possible incomplete row."""
    df = df.copy()
    df.drop_duplicates(inplace=True)
    mframe = df.groupby("frame").count()
    try:
        frame_from_rm = list(mframe[mframe["idx"] < mframe.loc[0, "idx"]].index)[0]
        df.drop(df[df["frame"] >= frame_from_rm].index, inplace=True)
    except IndexError:
        pass

    return df


class read_modes_vib:
    """Class defining the vibrational modes obtained in STAMP."""

    def __init__(self, file):
        """Initialize the class by reading the file."""
        self.file = file
        self._read_file()
        
    def _read_file(self):
        file = self.file
        tab = {}
        
        freq = re.compile(r"""
        ^nu\[(?P<idf>\d+)/\d+/\d+\]=(?P<f>\d+.\d+)
        """, re.X)
        
        coord = re.compile(r"""
        ^(?P<atsb>\w+)\s+
        (?P<x>[+-]?\d+\.\d+)\s+
        (?P<y>[+-]?\d+\.\d+)\s+
        (?P<z>[+-]?\d+\.\d+)\s+
        \w+\s+
        (?P<vx>[+-]?\d+\.\d+)\s+
        (?P<vy>[+-]?\d+\.\d+)\s+
        (?P<vz>[+-]?\d+\.\d+)\s+
        """, re.X)
        
        nfreq = 0
        with open(file, "r") as F:
            for line in F:
                array = line.split()
                if len(array) == 1:
                    nfreq += 1
            
                elif freq.match(line):
                    m = freq.match(line)
                    tab[int(m.groupdict()['idf'])] = dict()
                    tab[int(m.groupdict()['idf'])]['freq'] = float(m.groupdict()['f'])
                    tab[int(m.groupdict()['idf'])]["coord"] = []
            
                elif coord.match(line):
                    m = coord.match(line)
                    tab[nfreq]["coord"].append(m.groupdict())
                
        for i in tab:
            coord = pd.DataFrame(tab[i]["coord"])
            coord = coord.astype({
                'x': np.float64,
                'y': np.float64,
                'z': np.float64,
                'vx': np.float64,
                'vy': np.float64,
                'vz': np.float64})
            coord["atsb"] = coord["atsb"].apply(lambda x: x.upper())
            tab[i]["coord"] = coord
            
        self.table = tab

        modo = []
        freq = []
        for i in self.table:
            modo.append(i)
            freq.append(self.table[i]["freq"])
    
        self.modes_vib_freqs = pd.DataFrame({"n": modo, "freq": freq})

    def show_vectors(self, n, amplitud=2.):
        """Display the molecule with the vectors."""
        tab = self.table

        print("Frequency {} cm-1".format(tab[n]["freq"]))
        coords = tab[n]["coord"].loc[:, ["x", "y", "z"]].values

        atoms = list(tab[n]["coord"].loc[:, ["atsb"]].values.T[0])

        mol = Molecule(atoms, coords)

        import nglview as nv
        view = nv.show_pymatgen(mol)
        view.update_representation(opacity=0.4, camara='orthographic')
        for i in tab[n]["coord"].index:
            po = tab[n]["coord"].loc[i, ["x", "y", "z"]].values
            pf = po + amplitud*tab[n]["coord"].loc[i, ["vx", "vy", "vz"]].values
            view.shape.add_arrow(po, pf, [0, 0.4, 0.8], 0.2)
    
        return view


def read_geometry_file(file, freq, t0=0.0):
    """Read molecular geometry archives."""
    geom = pd.read_csv(
        file,
        index_col=0
    )
    geom["time"] = (geom.index.astype(float) - 1) * freq + t0
    
    return geom
