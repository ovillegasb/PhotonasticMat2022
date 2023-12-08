#!/usr/bin/env python
# -*- coding=utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import mdtraj as md
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput
import os
import glob

SOLVETS = {
    "THF": "THF",
    "PB": "CarbonTetraChloride"
}


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="geomSampling_2gaus",
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
        type=str,
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
        help="Output folder.",
        type=str,
        default="SamplesTDDFT"
    )

    analysis = parser.add_argument_group(
        "\033[1;36mAnalysis options\033[m")

    analysis.add_argument(
        "-b", "--b",
        help="Time of first frame to read from trajectory (default unit ps).",
        type=float,
        default=0
    )

    analysis.add_argument(
        "-e", "--e",
        help="Time of last frame to read from trajectory (default unit ps).",
        type=float,
        default=None
    )

    analysis.add_argument(
        "--init_t",
        help="Initial simulation time (default unit ps).",
        type=float,
        default=0.0
    )

    analysis.add_argument(
        "-dt", "--dt",
        help="Frame generation step (default unit ps).",
        type=float,
        default=1.0
    )

    analysis.add_argument(
        "-r", "--resid",
        help="resid number.",
        type=int,
        default=0
    )

    analysis.add_argument(
        "-sol", "--solvent",
        help="Additional comments for solvent, example: Tetrahidrofuran and CarbonTetraChloride",
        type=str,
        default=""
    )

    analysis.add_argument(
        "-isomer", "--isomer",
        help="Isomer type.",
        type=str,
        default=""
    )

    analysis.add_argument(
        "-N", "--NSamples",
        help="Number of samples to take.",
        type=int,
        default=100
    )

    return vars(parser.parse_args())


def clean_folder(folder):
    """Remove all files in a folder."""
    files = glob.glob(f"{folder}/*")
    for f in files:
        os.remove(f)


def read_xyzTraj(trajXYZ):
    """Return the information of an XYZ trajectory file."""
    xyz_list = []
    with open(trajXYZ, "r") as XYZ:
        nline_tot = 0
        nsaved = 0
        for line in XYZ:
            line = line.replace("\n", "").split()
            if len(line) == 1:
                nline_tot = int(line[0])
                nsaved = 0
                frame_values = []
                _ = XYZ.readline()
                continue

            elif len(line) == 4:
                frame_values.append(np.array(line))
                nsaved += 1

            if nsaved == nline_tot:
                frame = np.array(frame_values)
                atsb = frame[:, 0].astype(str)
                x = frame[:, 1].astype(np.float64)
                y = frame[:, 2].astype(np.float64)
                z = frame[:, 3].astype(np.float64)
                frame = pd.DataFrame({
                    "atsb": atsb,
                    "x": x,
                    "y": y,
                    "z": z
                })
                xyz_list.append(frame)

    return xyz_list


def main():
    """Run program."""
    args = options()

    out = args["out"]
    top = args["topol"]
    assert top is not None or args["traj"].endswith("xyz"), "You have to define a topology file."
    trj = args["traj"]
    assert trj is not None, "You have to define a trajectory file."
    isomer = args["isomer"]
    assert isomer is not None, "You have to define a isomer type."

    resid = args["resid"]
    out += "_r{}".format(resid)

    b = float(args["b"])
    if b != 0.0:
        out += "_{:.1f}_to".format(b)

    try:
        e = float(args["e"])
    except TypeError:
        out += ""
        e = -1
    else:
        out += "_{:.1f}".format(e)

    solvent = args["solvent"]
    if len(solvent) > 0:
        out += "_{}".format(solvent)

    if len(isomer) > 0:
        out += "_{}".format(isomer)

    samples_number = args["NSamples"]
    out += "_N{}".format(samples_number)

    if not os.path.exists(out):
        os.mkdir(out)
    else:
        clean_folder(out)

    print("File top:", top)
    print("File trj:", trj)

    # delta t regular
    dt = args["dt"]
    init_t = args["init_t"]

    print("Time >=", b, "ps", end=" - ")
    b -= init_t
    begin = int(b / dt)
    print("frame b:", begin)

    if e != -1:
        print("Time <", e, "ps", end=" - ")
        e -= init_t
        end = int(e / dt)
    else:
        print("frame e: Last")
        end = None
    print("Resid:", resid)

    #
    route_parameters = {}
    route_parameters["TD"] = "(NStates=6)"

    if len(solvent) != 0:
        try:
            route_parameters["SCRF"] = "(Solvent=%s)" % SOLVETS[solvent]
        except KeyError:
            print("ERROR solvent not recognized in the list")
            exit()

    frames_sample = []

    if top is not None and not trj.endswith("xyz"):
        assert top.endswith("gro"), "You have to define a topology gro file."
        # Reads the system trajectory
        t = md.load(trj, top=top)
        print("Total frames in trajectory:", t.n_frames)

        traj = t[begin:end]
        print("Frames selected:", traj.n_frames)
        frames_sample = np.random.choice(
            range(len(traj)), samples_number, replace=False
        )

        # Extracts the system topology
        top = traj.topology
        table, bonds = top.to_dataframe()

        # The indices of the atoms of interest are selected.
        iat_res = top.select("resid %s" % resid)

    elif trj.endswith("xyz"):
        # Read the system topolgy
        xyzTraj = read_xyzTraj(trj)
        print("Total frames in trajectory:", len(xyzTraj))

        traj = xyzTraj[begin:end]
        print("Frames selected:", len(traj))
        frames_sample = np.random.choice(
            range(len(traj)), samples_number, replace=False
        )

        atoms = traj[0]["atsb"].values
        xyz = []
        for frame in traj:
            xyz.append(frame.loc[:, ["x", "y", "z"]].values)
        xyz = np.array(xyz)

    # for frame to gaussian input
    for i in frames_sample:

        if not trj.endswith("xyz"):
            atoms = table.loc[iat_res, "element"].values
            coord = traj.xyz[:, iat_res][i] * 10.0  # to angstroms
        else:
            coord = xyz[i]
    
        mol = Molecule(atoms, coord)
        try:
            iso = isomer[0].upper() + "_"
        except IndexError:
            iso = ""

        if len(solvent) == 0:
            sol = ""
        else:
            sol = solvent + "_"
    
        # save file gaussian
        GaussianInput(
            mol,
            charge=0,
            spin_multiplicity=1,
            title="Azobenzene %s Sampling in %s - sample %d" % (isomer, solvent, i),
            functional="Cam-B3LYP",
            basis_set="6-311+g(d,p)",
            route_parameters=route_parameters,
            link0_parameters={
                "%mem": "8GB",
                "%nprocshared": "12"
            }
        ).write_file("%s/azo_%s%s%005d.com" % (out, iso, sol, i))

    print("Output folder: \033[1;36m{:>15}\033[m".format(out))


if __name__ == '__main__':
    main()
