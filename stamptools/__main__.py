#!/bin/env python

"""Using stamptools directly in terminal."""

import argparse
import os
from stamptools.stamp import STAMP
from stamptools.analysis import save_system, load_system
from stamptools.analysis import traj_analysis
from stamptools.analysis import traj_center_mass, get_distances_from
from stamptools.analysis import gen_centered_traj, mol_traj_analysis
from stamptools.analysis import mol_traj_cut_distance, get_angles_distance
from molcraft import clusters
from stamptools.stamptools import clean_data
import pandas as pd


TITLE = """\033[1;36m
   _____ _______       __  __ _____ _______ ____   ____  _       _____ 
  / ____|__   __|/\\   |  \\/  |  __ \\__   __/ __ \\ / __ \\| |     / ____|
 | (___    | |  /  \\  | \\  / | |__) | | | | |  | | |  | | |    | (___  
  \\___ \\   | | / /\\ \\ | |\\/| |  ___/  | | | |  | | |  | | |     \\___ \\ 
  ____) |  | |/ ____ \\| |  | | |      | | | |__| | |__| | |____ ____) |
 |_____/   |_/_/    \\_\\_|  |_|_|      |_|  \\____/ \\____/|______|_____/ 
\033[m
Module created to study systems generated by STAMP.

Author: Orlando VILLEGAS
Date: 2022-09-08
Stamp: v4.220721

Usage:

    Create object:
    python -m stamptools -d DONNEES.in

    Thermo:
    python -m stamptools -l -p T P Etot

    Save traj mol 0:
    python -m stamptools -l --mol 0

    Save poly information:
    python -um stamptools -l --poly > pol.log &

    Analysis of distances with respect to a resid:
    python -um stamptools -l --centerm -mref 0 > dist.log &

    System center from a reference:
    python -m stamptools -l --centered_traj --rcutoff 1.0 -mref 0

    Molecules around from a distance cut
    python -m stamptools -l --mol_d_aa -mref 0 --out_folder distances --rcuto\
ff 0.5

    On server:
    python -um stamptools -l > log &

    Combined analysis:
    python -m stamptools -l --centered_traj --rcutoff 1.0 -mref 0\
 --out_folder test0 --centerm

    python -m stamptools -l --poly --centerm -mref 0

    Angles distances:
    python -m stamptools -l --centerm --ref_plane -mref 0 -atref 11 12 13

"""


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="stamptools",
        usage="%(prog)s [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    fileinput = parser.add_argument_group(
        "\033[1;36mInitial settings\033[m")

    fileinput.add_argument(
        "-d", "--donnees",
        help="Specifies the system DONNEES file.",
        default=None
    )

    fileinput.add_argument(
        "--dataStamp",
        help="Specify data file, default is Stamp.dat.",
        default="Stamp.dat"
    )

    fileinput.add_argument(
        "-l", "--load",
        help="Loads a system status file.",
        action="store_true"
    )

    analysis = parser.add_argument_group(
        "\033[1;36mAnalysis options\033[m")

    analysis.add_argument(
        "-p", "--plots",
        help="""Saves graphs of thermodynamic parameters.

        Thermodynamic parameters: T, P or Etot, Several can be chosen
        at the same time.""",  # action="store_true"
        nargs="+",
        default=None
    )

    analysis.add_argument(
        "--mol",
        help="Analyze trajectory of a particular molecule, use resid.",
        type=int,
        default=None
    )

    analysis.add_argument(
        "--poly",
        help="Analyze the shape and size of the polymers present.",
        action="store_true"
    )

    analysis.add_argument(
        "--centerm",
        help="Analyze the of center of mass from molecules present.",
        action="store_true"
    )

    analysis.add_argument(
        "-mref", "--mref",
        help="Reference molecule.",
        type=int,
        default=None
    )

    analysis.add_argument(
        "-atref", "--atref",
        help="References atoms.",
        type=int,
        nargs="+",
        default=None
    )

    analysis.add_argument(
        "--centered_traj",
        help="Generates a trajectory using a reference molecule.",
        action="store_true"
    )

    analysis.add_argument(
        "--rcutoff",
        help="Cut-off distance with respect to the reference molecule.",
        type=float,
        default=1.0
    )

    analysis.add_argument(
        "--out_folder",
        help="Name of the output directory.",
        type=str,
        default="centered_traj"
    )

    analysis.add_argument(
        "--mol_d_aa",
        help="Extract the structure of mol around a reference mol using atom-a\
tom distance.",
        action="store_true"
    )

    analysis.add_argument(
        "--mol_d",
        help="Analyze center of mass distance from a referece, use resid.",
        action="store_true"
    )

    analysis.add_argument(
        "--reset",
        help="Re-collect data.",
        action="store_true"
    )

    analysis.add_argument(
        "--ref_plane",
        help="The angle between the normal vector to a plane formed by three a\
toms centered at the center of mass of the reference molecule and the center o\
f mass is analyzed.",
        action="store_true"
    )

    return vars(parser.parse_args())


print(TITLE)
args = options()
# print(args)

if not args["load"] and args["donnees"]:
    system = STAMP(donnees=args["donnees"], data=args["dataStamp"])
    save_system(system)

elif args["load"]:
    print("The system status will be loaded")
    system = load_system("system.chk")

else:
    print("The state of the system must be defined.")
    exit()


# Others options:

if args["plots"]:
    system.save_plots(args["plots"])

if isinstance(args["mol"], int):
    resid = args["mol"]
    print("Resid:", resid)
    mol_ndx = system.atoms_per_mol[resid]
    mol_traj_analysis(
        resid,
        mol_ndx,
        system.connectivity,
        system.traj,
        system.box
    )

if args["poly"]:
    # begin frame
    b = 0
    if args["reset"]:
        os.remove("polymers.csv")
        print("File polymers.csv removed.")

    if os.path.exists("polymers.csv"):
        dat = pd.read_csv("polymers.csv")
        dat = clean_data(dat)
        dat.to_csv("polymers.csv", index=False)
        frames_readed = list(pd.unique(dat["frame"]))
        # update file list
        system.update_xyz()

        if len(frames_readed) == len(system.XYZs):
            print("The number of XYZ files is equal to the number of "\
                  "frames analyzed.")
            b = len(frames_readed)

        elif len(frames_readed) < len(system.XYZs):
            print("The number of XYZ files is greater than the number of "\
                  "files analyzed.")
            b = len(frames_readed)
            save_system(system)

    if b == 0 and not args["reset"]:
        args["reset"] = True

    traj_analysis(
        system.atoms_per_mol,
        system.topology, 
        system.traj,
        system.box,
        system.connectivity,
        clusters.GyrationTensor,
        b,
        args["reset"]
    )
    # save information in file
    print("file polymers.csv saved.")

if args["centerm"]:
    # begin frame
    b = 0
    if args["reset"]:
        os.remove("mol_cmass.csv")
        print("File mol_cmass.csv removed.")

    if os.path.exists("mol_cmass.csv"):
        dat = pd.read_csv("mol_cmass.csv")
        dat = clean_data(dat)
        dat.to_csv("mol_cmass.csv", index=False)
        frames_readed = list(pd.unique(dat["frame"]))
        # update file list
        system.update_xyz()

        if len(frames_readed) == len(system.XYZs):
            print("The number of XYZ files is equal to the number of "\
                  "frames analyzed.")
            b = len(frames_readed)

        elif len(frames_readed) < len(system.XYZs):
            print("The number of XYZ files is greater than the number of "\
                  "files analyzed.")
            b = len(frames_readed)
            save_system(system)

    if b == 0 and not args["reset"]:
        args["reset"] = True

    traj_center_mass(
        system.traj,
        system.atoms_per_mol,
        system.topology,
        system.box,
        system.connectivity,
        b,
        args["reset"]
    )
    # save information in file
    print("file mol_cmass.csv saved.")

if args["mol_d_aa"]:
    mol_traj_cut_distance(
        system,
        ref=args["mref"],
        rcutoff=args["rcutoff"],
        out_folder=args["out_folder"]
    )

if args["mol_d"]:
    print("Resid:", args["mref"])
    get_distances_from(args["mref"], system.box)
    print("file mol_dist_from_{}.csv saved.".format(args["mref"]))

if args["centered_traj"]:
    # distances file
    mol_dist = pd.read_csv("mol_dist_from_{}.csv".format(args["mref"]))
    mol_dist["distance"] = mol_dist["distance"] * 0.1  # to nm

    # load center of mass file
    c_mass = pd.read_csv("mol_cmass.csv")

    gen_centered_traj(
        system,
        mol_dist,
        c_mass,
        rcutoff=args["rcutoff"],
        ref=args["mref"],
        out_folder=args["out_folder"])

if args["ref_plane"]:
    print("Analysis using a plane")

    # Reference molecule using mref
    mref = args["mref"]
    print("Resid:", mref)

    # References atoms indexs
    atref = args["atref"]
    print("Index:", " ".join([str(i) for i in atref]))

    # FIle traj center of mass
    file = "mol_cmass.csv"

    # Coordinates
    traj = system.traj
    box = system.box

    get_angles_distance(mref, atref, box, traj, file)
    print(f"file mol_angles_d_{mref}.csv saved.")
