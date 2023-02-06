#!/bin/env python

"""Module to define the STAMP class."""

import glob
import time
import os
import numpy as np
from pathlib import Path
from scipy.constants import N_A
from stamptools.stamptools import read_donnees
from stamptools.analysis import load_data, read_fatomes, save_plot, load_log
from molcraft import structure
from multiprocessing import Pool

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


class STAMP:
    """Object describing a system and its features."""

    def __init__(self, donnees, data="Stamp.dat", ensamble="LNVT"):
        """Initialize the class when loading the calculation information."""
        # Home path
        hw_path = os.path.split(os.path.abspath(donnees))[0]
        self.hw_path = hw_path

        # Load donnees informations
        self.donnees = read_donnees(donnees)

        # Define Ensemble
        self.ensemble = self.donnees["Ensemble"]

        # File fatomes
        self.fatomes = self.donnees["FichierAtomes"]

        # Load results from Stamp.dat
        self.data = load_data(os.path.join(hw_path, data), t=ensamble)

        # Load log information
        try:
            self.time_per_frame = load_log(os.path.join(hw_path, "Stamp.log"))
        except FileNotFoundError:
            self.time_per_frame = None

        # Reading FAtome
        topology, box, connects = read_fatomes(
            os.path.join(hw_path, self.fatomes)
        )

        self.topology = topology
        self.box = box
        self.connects = connects

        # searching conectivity
        self.connectivity = self._get_connectivity()

        # Files xyz
        self.XYZs = self._xyz_list()
        # if load_traj:
        #     self._load_traj()

    def save_plots(self, args):
        """Save plot for parameters."""
        dat = self.data
        for a in args:
            if a in setplots:
                save_plot(
                    dat["I"], dat[a],
                    name=setplots[a]["name"],
                    color=setplots[a]["color"],
                    xlb=setplots[a]["xlb"],
                    ylb=setplots[a]["ylb"]
                )

                print("plots {}.png saved".format(setplots[a]["name"]))

            else:
                print(f"Option {a} has not been configured.")

    @property
    def box_in_time(self):
        """Box [ang, ang, ang] in times from xyz files."""
        boxs = []
        if self.time_per_frame is not None:
            for j in self.time_per_frame.index:
                with open(self.XYZs[j], "r") as xyz:
                    for i, line in enumerate(xyz):
                        if i == 1:
                            line = line.replace("\n", "").split()
                            boxs.append(line)
                            break
        else:
            for i in self.XYZs:
                with open(i, "r") as xyz:
                    for i, line in enumerate(xyz):
                        if i == 1:
                            line = line.replace("\n", "").split()
                            boxs.append(line)
                            break

        return np.array(boxs).astype(np.float)


    @property
    def vol(self):
        """Volume (nm3) from xyz file stamp."""
        boxs = self.box_in_time
        return boxs[:, 0] * boxs[:, 1] * boxs[:, 2] * 1e-3

    @property
    def dens(self):
        """Density from xyz files stamp."""
        massT = self.topology.mass.sum() / N_A

        return massT / self.vol / 1e-27

    @property
    def dens_in_time(self):
        dens = self.time_per_frame
        dens["Dens"] = self.dens
        return dens

    def _xyz_list(self):
        """Load list of file xyz."""
        return sorted(glob.glob(f"{self.hw_path}/XYZ/PasDeCalcul_Iteration*.xyz"))

    def update_xyz(self):
        """Reload list of file xyz and traj."""
        self.XYZs = self._xyz_list()
        # self._load_traj()

    def _load_traj(self, b=0, e=None):
        """Load trajectory system."""
        traj = list()
        t0 = time.time()
        print("Loading the system trajectory", end=" - ")

        with Pool() as pool:
            for xyz in pool.map(structure.load_xyz, self.XYZs[b:e]):
                traj.append(xyz)

        self._traj = traj
        tf = time.time()
        print(f"Number of frames: {len(traj)}", end=" - ")
        print(f"done in {tf-t0:.2f} s")

    def get_traj(self, b=0, e=None):
        """Trajectory of system."""
        self._load_traj(b, e)

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
