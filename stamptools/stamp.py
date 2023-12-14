#!/bin/env python

"""Module to define the STAMP class."""

import glob
import time
import os
import numpy as np
from scipy.constants import N_A
from stamptools.stamptools import read_donnees
from stamptools.analysis import load_data, read_fatomes, load_log
from molcraft import structure
from multiprocessing import Pool
from stamptools.gmxtools import read_xtc

NUMPROC = 12

setplots = {
    "T": {
        "xlb": "time (ps)", "ylb": "Temperature (K)",
        "name": "temp", "color": "r"
        },
    "P": {
        "xlb": "time (ps)", "ylb": "Pressure (bar)",
        "name": "press", "color": "green"
        },
    "Etot": {
        "xlb": "time (ps)", "ylb": "Total energy (kj/mol)",
        "name": "etot", "color": "b"
        }
}

toElements = {
    "OHopls": "oh",
    "HOopls": "ho",
    "CTO": "CT",
    "H1O": "HT"
}


def change_toType(at):
    try:
        return toElements[at]
    except KeyError:
        return at


class STAMP:
    """Object describing a system and its features."""

    def __init__(self, donnees, fatomes=None, data="Stamp.dat", ensamble="LNVT", loadConnect=True, traj_type="XYZ", use_xyz=True, nproc=NUMPROC):
        """Initialize the class when loading the calculation information."""
        # Home path
        hw_path = os.path.split(os.path.abspath(donnees))[0]
        self.hw_path = hw_path

        # Number of cpus
        self.NUMPROC = nproc

        # Load donnees informations
        self.donnees = read_donnees(donnees)

        # Define Ensemble
        self.ensemble = self.donnees["Ensemble"]

        # File fatomes
        if fatomes is None:
            self.fatomes = self.donnees["FichierAtomes"]
        else:
            self.fatomes = fatomes

        # Load results from Stamp.dat
        self.data = load_data(os.path.join(hw_path, data), t=ensamble)

        # Load log information
        try:
            self.time_per_frame, self.status = load_log(os.path.join(hw_path, "Stamp.log"), use_xyz=use_xyz, traj_type=traj_type)
        except FileNotFoundError:
            self.time_per_frame = None

        # Reading FAtome
        if loadConnect:
            topology, box, connects = read_fatomes(
                os.path.join(hw_path, self.fatomes)
            )

            topology["atsb"] = topology["atsb"].apply(change_toType)

            self.topology = topology
            self.box = box
            self.connects = connects

            # searching conectivity
            self.connectivity = self._get_connectivity()

        if traj_type == "XYZ":
            # Files xyz
            self.XYZs = self._xyz_list()
            print("Number of frames", len(self.XYZs))
        elif traj_type == "GRO":
            # Files gro
            self.GROs = self._gro_list()
            print("Number of frames", len(self.GROs))
        elif traj_type == "XTC":
            self.top = os.path.join(hw_path, "XTC/confout.gro")
            self.xtc = os.path.join(hw_path, "XTC/traj_comp.xtc")
            self._box_xtc = None
        self.b_frame = 0
        self.e_frame = None
        self.i_frame = 1
        self.traj_type = traj_type
        # if load_traj:
        #     self._load_traj()

    @property
    def Natoms(self):
        """Total number of atoms in the system."""
        return len(self.topology)

    @property
    def box_in_frame(self):
        """Box [ang, ang, ang] in times from xyz files."""
        boxs = []
        if self.traj_type == "XYZ":
            e = self.e_frame if self.e_frame is not None else len(self.XYZs)
            for file in self.XYZs[self.b_frame:e+1:self.i_frame]:
                with open(file, "r") as xyz:
                    for i, line in enumerate(xyz):
                        if i == 1:
                            line = line.replace("\n", "").split()
                            boxs.append(line)
                            break
            return np.array(boxs).astype(np.float64)

        elif self.traj_type == "GRO":
            e = self.e_frame if self.e_frame is not None else len(self.GROs)
            for file in self.GROs[self.b_frame:e+1:self.i_frame]:
                with open(file, 'r') as f:
                    last_line = f.readlines()[-1]
                    last_line = last_line.replace("\n", "").split()
                    boxs.append(last_line)
            return np.array(boxs).astype(np.float64) * 10.

        elif self.traj_type == "XTC":
            assert self._box_xtc is not None, "To read the box information using XTC you must load the system trajectory. STAMP.get_traj"
            boxs = self._box_xtc
            return np.array(boxs).astype(np.float64) * 10.

    @property
    def vol(self):
        """Volume (nm3) from xyz file stamp."""
        boxs = self.box_in_frame
        return boxs[:, 0] * boxs[:, 1] * boxs[:, 2] * 1e-3

    @property
    def dens(self):
        """Density from xyz files stamp."""
        massT = self.topology.mass.sum() / N_A

        return massT / self.vol / 1e-27

    @property
    def dens_in_time(self):
        """Density in time from boxs."""
        dens = self.time_per_frame
        dens["Dens"] = self.dens
        return dens

    def _xyz_list(self):
        """Load list of file xyz."""
        return sorted(glob.glob(f"{self.hw_path}/XYZ/PasDeCalcul_Iteration*.xyz"))

    def _gro_list(self):
        """Load list of file gro."""
        return sorted(glob.glob(f"{self.hw_path}/GRO/PasDeCalcul_*.gro"))

    def update_xyz(self):
        """Reload list of files and traj."""
        if self.traj_type == "XYZ":
            self.XYZs = self._xyz_list()
        elif self.traj_type == "GRO":
            self.GROs = self._gro_list()

    def _load_traj(self, b=None, e=None, i=None):
        """Load trajectory system."""
        traj = list()
        t0 = time.time()
        print("Loading the system trajectory", end=" - ")

        if self.traj_type == "XYZ":
            XYZs = self.XYZs
            e = e if e is not None else len(XYZs)
            with Pool(processes=self.NUMPROC) as pool:
                for xyz in pool.map(structure.load_xyz, XYZs[b:e+1:i]):
                    traj.append(xyz)
        elif self.traj_type == "GRO":
            GROs = self.GROs
            e = e if e is not None else len(GROs)
            with Pool(processes=self.NUMPROC) as pool:
                for gro in pool.map(structure.load_gro, GROs[b:e+1:i]):
                    traj.append(gro)
        elif self.traj_type == "XTC":
            b = b if b is not None else 0
            i = i if i is not None else 1
            params = {
                "top": self.top,
                "xtc": self.xtc,
                "interval": i,
                "b": b,
                "e": e
            }
            traj, boxs = read_xtc(**params)
            self._box_xtc = boxs

        self._traj = traj
        self.b_frame = b
        self.e_frame = e
        self.i_frame = i
        tf = time.time()
        print(f"Number of frames: {len(traj)}", end=" - ")
        print(f"done in {tf-t0:.2f} s")

    def get_traj(self, b=None, e=None, i=None):
        """Trajectory of system."""
        self._load_traj(b, e, i)
        return self._traj        

    def _get_connectivity(self):
        """Brings system connectivity."""
        t0 = time.time()
        print("Searching connectivity", end=" - ")
        conn = structure.connectivity()
        conn.define_atoms(self.topology)
        conn.read_dict(self.connects)

        tf = time.time()
        print(f"done in {tf-t0:.2f} s")

        return conn

    @property
    def atoms_per_mol(self):
        """List atoms per molecule."""
        return self.connectivity.atomsMOL

    def monomers_in_chains(self, natoms_in_monomer=10):
        """Back a list of monomers in polymer chains."""
        monomers = []
        atoms_per_mol = self.atoms_per_mol
        for imol in atoms_per_mol:
            natoms = atoms_per_mol[imol]["Natoms"]
            natoms_in_monomer = 10
            if natoms > 1:
                monomers.append(natoms / natoms_in_monomer)

        return monomers
