#!/bin/env python

"""RDF analysis all atoms between systems."""

import argparse
import glob
import numpy as np
from stamptools.analysis import read_fatomes, PBC_distance, progress
from molcraft.structure import connectivity, load_xyz
import sys
from scipy.spatial.distance import cdist
import time
import pandas as pd
import matplotlib.pyplot as plt

"""
# Folder and systems
systems = ["0", "1", "2", "3", "4", "0_long"]
isomer = "azoC_procedure"
home = "/home/ovillegas"
fatomes = f"{home}/{isomer}/6_prod_0/FAtomes.inp_6_prod"
name = "RDF_all_cis"

# RDFs options
resid = 0
# Define the bins
rmin = 0.0
rmax = 3.0
binwidth = 0.05
bins = np.arange(rmin, rmax + binwidth, binwidth)

# Reading FAtome
topology, box, connects = read_fatomes(fatomes)
vol = box[0] * box[1] * box[2] * 0.1**3

conn = connectivity()
conn.define_atoms(topology)
conn.read_dict(connects)
atoms_per_mol = conn.atomsMOL

atoms_ref = atoms_per_mol[resid]["index"]
atoms_env = []
for i in atoms_per_mol:
    if i != resid:
        atoms_env += atoms_per_mol[i]["index"]

print("resid:      ", resid)
print("N atoms ref:", len(atoms_ref))
print("N atoms env:", len(atoms_env))


XYZs = []
for s in systems:
    XYZs += glob.glob(f"{home}/{isomer}/6_prod_{s}/XYZ/*.xyz")

print("N of frames:", len(XYZs))

frame_distances = {}
print(
    "Memory in frame_distances before traj analysis:",
    sys.getsizeof(frame_distances),
    "bytes"
)

start = time.perf_counter()

for i, frame in enumerate(XYZs):
    coord = load_xyz(frame, warning=False)
    # print(frame)
    atoms_xyz_ref = coord.loc[atoms_ref, ["x", "y", "z"]].values
    atoms_xyz_env = coord.loc[atoms_env, ["x", "y", "z"]].values
    # print(atoms_xyz_ref.shape)
    # print(atoms_xyz_env.shape)
    frame_distances[frame] = cdist(
            atoms_xyz_ref,
            atoms_xyz_env,
            lambda a, b: PBC_distance(a, b, box)
        )
    # if i == 50:
    #     break

end = time.perf_counter()

print("frame_distances computed: " + str(round(end - start, 2)) + "s")
print("Nomber of frames added:", len(frame_distances))
print(
    "Memory in frame_distances:",
    sys.getsizeof(frame_distances),
    "bytes"
)


"""

TITLE = """
Module created to calculate the RDF of a molecular system.

Author: Orlando VILLEGAS
Date: 2023-01-23

Usage:

    python compute_RDF.py -s FAtomes.inp_6_prod -ref resid 0 -sel resid 0 -f \
XYZ/*.xyz --show_plot -rmax 1. -bin 0.01

"""


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="compute_RDF",
        usage="%(prog)s [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    fileinput = parser.add_argument_group(
        "\033[1;36mInitial settings\033[m"
    )

    fileinput.add_argument(
        "-f", "--traj",
        help="Input trajectory or single configuration: xtc, gro, *.xyz.",
        nargs="+",
        default=None
    )

    fileinput.add_argument(
        "-s", "--topol",
        help="Input structure, system topology: gro, fatomes.",
        type=str,
        default=None
    )

    fileinput.add_argument(
        "-o", "--out",
        help="Computed RDFs (rdf.csv).",
        type=str,
        default="rdf.csv"
    )

    analysis = parser.add_argument_group(
        "\033[1;36mAnalysis options\033[m")

    analysis.add_argument(
        "-bin", "--bin",
        help="Bin width (nm).",
        type=float,
        default=0.05
    )

    analysis.add_argument(
        "-rmax", "--rmax",
        help="Largest distance (nm) to calculate.",
        type=float,
        default=3.0
    )

    analysis.add_argument(
        "-rmin", "--rmin",
        help="Shortest distance (nm) to calculate.",
        type=float,
        default=0.0
    )

    analysis.add_argument(
        "--show_plot",
        help="Plots the RDF obtained.",
        action="store_true",
        default=False
    )

    selections = parser.add_argument_group(
        """\033[1;36mSelection options\033[m
The selections work with the keyword resid, followed by the number or \
numbers of residues to be studied. You can also use conditions, exampl\
es:
resid N1 ... Nn
resid not N1 ... Nn
resid all
""")

    selections.add_argument(
        "-ref", "--ref",
        help="Reference selection for RDF computation.",
        nargs="+",
        type=str
    )

    selections.add_argument(
        "-sel", "--sel",
        help="Selections to compute RDFs for from the reference. ",
        nargs="+",
        type=str
    )

    return vars(parser.parse_args())


def main():
    """Execute the program."""
    print(TITLE)
    args = options()
    # print(args)

    # RDFs options
    rmin = args["rmin"]
    rmax = args["rmax"]
    binwidth = args["bin"]
    print("Start Program:")
    print(f"    rmin {rmin} nm")
    print(f"    rmax {rmax} nm")
    print(f"    bin {binwidth} nm")
    bins = np.arange(rmin, rmax + binwidth, binwidth)
    # print(bins)

    # Selection: ["resid", "N1", ..., "Nn"]
    ref = args["ref"]
    sel = args["sel"]
    assert ref is not None and sel is not None, "You must define a selection."
    assert ref[0] == "resid" and sel[0] == "resid", "You must add a selection \
keyword."

    # Trajectory
    traj = args["traj"]
    assert traj is not None, "No trajectory file has been selected."

    # Frame distances structure
    frame_distances = {}

    topol = args["topol"]
    if topol is not None:
        # print(topol)
        # Reading FAtomes
        topology, box, connects = read_fatomes(topol)
        vol = box[0] * box[1] * box[2] * 0.1**3

        conn = connectivity()
        conn.define_atoms(topology)
        conn.read_dict(connects)
        # dict = {resid : {"Natoms": int, "index": [...]}}
        atoms_per_mol = conn.atomsMOL

        # references atomes
        atoms_ref = []
        if "resid" == ref[0]:
            for i in ref[1:]:
                atoms_ref += atoms_per_mol[int(i)]["index"]

        atoms_sel = []
        if "resid" == sel[0]:
            for i in sel[1:]:
                atoms_sel += atoms_per_mol[int(i)]["index"]
        
        print("N atoms ref:", len(atoms_ref))
        print("N atoms sel:", len(atoms_sel))
        # Number of frames read
        Nframes = len(traj)
        print(f"Number of frames: {Nframes}")
        start = time.perf_counter()

        for i, frame in enumerate(traj):
            porcent = i * 100 / Nframes
            print(f"{porcent:6.2f} % |{progress(porcent)}|")
            coord = load_xyz(frame, warning=False)
            atoms_xyz_ref = coord.loc[atoms_ref, ["x", "y", "z"]].values
            atoms_xyz_sel = coord.loc[atoms_sel, ["x", "y", "z"]].values

            frame_distances[frame] = cdist(
                atoms_xyz_ref,
                atoms_xyz_sel,
                lambda a, b: PBC_distance(a, b, box)
            )

        print(f"{100:6.2f} % |{progress(100)}|")

        end = time.perf_counter()

        print("frame_distances computed: " + str(round(end - start, 2)) + "s")
        print("Nomber of frames added:", len(frame_distances))
        print(
            "Memory in frame_distances:",
            sys.getsizeof(frame_distances),
            "bytes"
        )

    # RDFs
    # print(frame_distances)
    rdf_start = time.perf_counter()
    n_centers_ref = len(atoms_ref)
    n_centers_sel = len(atoms_sel)

    vol_per_sphere = vol / n_centers_sel
    vshell = 4 * np.pi * ((binwidth + bins)**3 - bins**3) / 3

    g_r = np.zeros(len(bins))
    for frame in frame_distances:
        for i, atom in enumerate(frame_distances[frame]):
            indexs = np.int64(atom * 0.1 / binwidth)
            indexs = indexs[indexs < len(bins)]
            for n in indexs:
                g_r[n] += 1

    n_frames = len(frame_distances.keys())
    g_r_norm = g_r * vol_per_sphere / vshell / n_frames / n_centers_ref

    # return g_r_norm, bins
    RDF = pd.DataFrame({
        "g_r": g_r_norm,
        "r": bins
    })

    rdf_end = time.perf_counter()
    print("RDF computed: " + str(round(rdf_end - rdf_start, 2)) + "s")

    file = args["out"]
    RDF.to_csv(file, float_format="%.6f", index=None)
    print(f"file {file} saved.")
    print("")

    if args["show_plot"]:
        fig, ax = plt.subplots()
        ax.plot(RDF["r"], RDF["g_r"])
        ax.set_xlim(rmin, rmax)
        ax.axhline(y=1, ls="--", color="gray")
        ax.set_xlabel("r (nm)")
        ax.set_ylabel("g(r)")
        plt.show()


if __name__ == '__main__':
    main()
