#!/bin/env python

"""
Script used to convert STAMP xyz trajectory to GROMACS gro trajectory.

Inputs options:

"""

import argparse
import os
import time
import glob
import numpy as np
from multiprocessing import Pool
from stamptools.stamp import STAMP
from stamptools.analysis import center_of_mass, minImagenC


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
            prog="xyz2gro",
            usage="%(prog)s [-options]",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="Enjoy the program!"
        )

    systemInput = parser.add_argument_group(
        "\033[1;36mInitial settings of system\033[m")

    systemInput.add_argument(
        "-d", "--donnees",
        help="Specifies the system DONNEES file.",
        default=None
    )

    return vars(parser.parse_args())


def save_gro(table, name, box, title="GRO FILE", time=0.0, out="."):
    """Save coordinate to file *.gro."""
    nat = len(table)
    gro = name

    GRO = open(f"{out}/{gro}.gro", "w", encoding="utf-8")
    GRO.write("{}, t= {:.3f}\n".format(title, time))
    GRO.write("%5d\n" % nat)

    for i in table.index:
        GRO.write("{:>8}{:>7}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
            str(table.loc[i, "resid"]) + table.loc[i, "resname"],
            table.loc[i, "atsb"],
            i + 1,
            table.loc[i, "x"] * 0.1 + box[0] / 20,
            table.loc[i, "y"] * 0.1 + box[1] / 20,
            table.loc[i, "z"] * 0.1 + box[2] / 20)
        )

    GRO.write("   {:.5f}   {:.5f}   {:.5f}\n".format(
        box[0] / 10,
        box[1] / 10,
        box[2] / 10))

    GRO.close()
    # print("\nSaved gro file: \033[1;36m%s\033[m writed\n" % gro)


def save_frame(i, frame, RESinfo, box, time, title, out):
    """Save information per a frame i."""
    frame["resid"] = frame.index.map(lambda x: RESinfo[x]["resid"])
    frame["resname"] = frame.index.map(lambda x: RESinfo[x]["resname"])
    save_gro(
        frame,
        "frame_%005d" % i,
        box,
        title,
        time=time,
        out=out)


def traj_to_xtc(traj, boxs, time_per_frame, RESinfo, out):
    """Convert xyz trajectory to xtc."""
    title = "PC in polymer"

    arguments = []
    for i, frame in enumerate(traj):
        arguments.append(
            (i, frame, RESinfo, boxs[i],
                time_per_frame.loc[i, "time"], title, out))

    with Pool() as pool:
        for _ in pool.starmap(save_frame, arguments):
            pass


def get_nopbc_mol(mol_ndx, confout, connectivity, box):
    """Return the structure complet fo a molecule."""
    mol_xyz = confout.loc[mol_ndx, :]
    mol_conn = connectivity.sub_connect(mol_ndx)
    mol_conn.update_coordinates(mol_xyz)
    mol_conn.noPBC(box, center=np.zeros(3))
    mol_conn = mol_conn.reset_nodes()
    mol_conn.simple_at_symbols(add_mass=True)
    return mol_conn.get_df(), mol_ndx


def clean_folder(folder):
    """Remove all files in a folder."""
    files = glob.glob(f"{folder}/*")
    for f in files:
        os.remove(f)


def main():
    """Run program."""
    args = options()
    out = "GRO_files"

    if not os.path.exists(out):
        os.mkdir(out)
    else:
        clean_folder(out)
    
    # Load System
    system = STAMP(donnees=args["donnees"])
    # ######## Function for correction of box
    # veloc_x = 0
    # try:
    #     veloc_x = float(system.donnees["PistonDvitesse"]) * 1e10 / 1e12 
    # except KeyError:
    #     pass
    # ########

    resid_PC = 0
    atoms_per_mol = system.atoms_per_mol
    for resid in atoms_per_mol:
        if resid == resid_PC:
            atoms_per_mol[resid]["resname"] = "PHO"
        elif len(atoms_per_mol[resid]["index"]) == 1:
            atoms_per_mol[resid]["resname"] = "PHO"
        else:
            atoms_per_mol[resid]["resname"] = "POL"

    RESinfo = {}
    for resid in atoms_per_mol:
        resname = atoms_per_mol[resid]["resname"]
        index = atoms_per_mol[resid]["index"]
        for i in index:
            RESinfo[i] = {"resid": resid, "resname": resname}

    time_per_frame = system.time_per_frame

    # get trajectory
    traj = system.get_traj()
    # # Is piston
    # if veloc_x != 0.0:
    #     for i, frame in enumerate(traj):
    #         frame["x"] = frame["x"] + time_per_frame.loc[i, "time"] * veloc_x

    # BOXs
    boxs = []
    for file in system.XYZs:
        with open(file, "r") as xyz:
            _ = xyz.readline()
            boxs.append(xyz.readline().split()[0:3])

    boxs = np.array(boxs).astype(np.float64)

    print("Saving GRO files", end=" - ")
    t0 = time.perf_counter()
    # convert to xtc
    traj_to_xtc(traj, boxs, time_per_frame, RESinfo, out)
    t1 = time.perf_counter()
    print("done in", str(round(t1 - t0, 2)) + "s")

    # Save a conf final with all atoms complets
    print("Generating a gro file without PBC", end=" - ")
    t0 = time.perf_counter()
    confout = traj[-1]
    confout["resid"] = confout.index.map(lambda x: RESinfo[x]["resid"])
    confout["resname"] = confout.index.map(lambda x: RESinfo[x]["resname"])

    box = boxs[-1]
    connectivity = system.connectivity
    pc_xyz = confout[confout["resname"] == "PHO"]
    pc_conn = connectivity.sub_connect(atoms_per_mol[resid_PC]["index"])
    # update coordinates
    pc_conn.update_coordinates(pc_xyz)
    # remove PBC
    pc_conn.noPBC(box)
    # Reset index and symbols, and add mass
    pc_conn = pc_conn.reset_nodes()
    pc_conn.simple_at_symbols(add_mass=True)

    new_pc_xyz = pc_conn.get_df()
    cm = center_of_mass(
            new_pc_xyz.loc[:, ["x", "y", "z"]].values,
            new_pc_xyz.loc[:, "mass"].values
        )

    new_pc_xyz["x"] -= cm[0]
    new_pc_xyz["y"] -= cm[1]
    new_pc_xyz["z"] -= cm[2]

    confout["x"] = confout["x"].apply(lambda x: minImagenC(cm[0], x, box[0]))
    confout["y"] = confout["y"].apply(lambda x: minImagenC(cm[1], x, box[1]))
    confout["z"] = confout["z"].apply(lambda x: minImagenC(cm[2], x, box[2]))

    arguments = []
    for i in atoms_per_mol:
        arguments.append(
            (atoms_per_mol[i]["index"], confout, connectivity, box)
        )

    with Pool() as pool:
        for new_mol_xyz, mol_ndx in pool.starmap(get_nopbc_mol, arguments):
            confout.loc[mol_ndx, "x"] = new_mol_xyz["x"].values
            confout.loc[mol_ndx, "y"] = new_mol_xyz["y"].values
            confout.loc[mol_ndx, "z"] = new_mol_xyz["z"].values

    save_gro(confout, "confout", box, out=out)
    t1 = time.perf_counter()
    print("done in", str(round(t1 - t0, 2)) + "s")

    os.chdir(out)
    os.system("cat frame_* > traj.gro")
    os.system("rm frame_*")
    os.system("echo 0 | gmx trjconv -f traj.gro -s confout.gro -o traj_comp.xtc")
    os.system("rm traj.gro")
    os.system("cd ..")
    print("Saved trajectory gro in traj_comp.xtc")


if __name__ == '__main__':
    main()
