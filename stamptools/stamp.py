#!/bin/env python

"""Module to define the STAMP class."""

import glob
import numpy as np
from scipy.constants import N_A
from stamptools.stamptools import read_donnees
from stamptools.analysis import load_data, read_fatomes, save_plot, traj_analysis
from molcraft import structure, clusters

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

    def __init__(self, donnees, data="Stamp.dat"):
        """Initialize the class when loading the calculation information."""
        self.donnees = read_donnees(donnees)
        # print(self.donnees)

        # Define Ensemble
        self.ensemble = self.donnees["Ensemble"]

        # File fatomes
        self.fatomes = self.donnees["FichierAtomes"]

        # Load results from Stamp.dat
        self.data = load_data(data)

        # Reading FAtome
        topology, box, connects = read_fatomes(self.fatomes)
        self.topology = topology
        self.box = box
        self.connects = connects

        self.connectivity = self._get_connectivity()

        # print(self.topology)
        # print(self.box)
        # print(self.connects)

        # Files xyz
        self.XYZs = sorted(glob.glob("./XYZ/PasDeCalcul_Iteration*.xyz"))

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
    def dens(self):
        """Density from xyz files stamp."""
        info = []
        for i in self.XYZs:
            with open(i, "r") as xyz:
                for i, line in enumerate(xyz):
                    if i == 1:
                        line = line.replace("\n", "").split()
                        info.append(line)
                        break

        info = np.array(info).astype(np.float)
        vol = info[:, 0] * info[:, 1] * info[:, 2] * 1e-27
        massT = self.topology.mass.sum() / N_A

        return massT / vol

    @property
    def traj(self):
        """Trajectory of system."""
        traj = list()

        for file in self.XYZs:
            xyz = structure.load_xyz(file, warning=False)
            traj.append(xyz)

        return traj

    def _get_connectivity(self):
        """Brings system connectivity."""
        conn = structure.connectivity()
        conn.define_atoms(self.topology)
        conn.read_dict(self.connects)

        return conn

    @property
    def atoms_per_mol(self):
        """List atoms per molecule."""
        return self.connectivity.atomsMOL

    def get_poly_info(self):
        """Save file information for size mol."""
        dfpoly = traj_analysis(
            self.atoms_per_mol,
            self.topology, 
            self.traj,
            self.box,
            self.connectivity,
            clusters.GyrationTensor
        )

        # save information in file
        dfpoly.to_csv("datapolyP.csv")
        print("file datapolyP saved.")
