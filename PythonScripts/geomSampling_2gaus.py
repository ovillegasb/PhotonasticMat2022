#!/usr/bin/env python
# -*- coding=utf-8 -*-

import argparse
import numpy as np
import mdtraj as md
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput
import os
import glob


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
        default="SamplingUV-Vis"
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
        "-r", "--resid",
        help="resid number.",
        type=int,
        default=0
    )

    analysis.add_argument(
        "-sol", "--solvent",
        help="Additional comments for solvent.",
        type=str,
        default=""
    )

    analysis.add_argument(
        "-isomer", "--isomer",
        help="Isomer type.",
        type=str,
        default=None
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


def main():
    """Run program."""
    args = options()

    out = args["out"]

    top = args["topol"]
    assert top is not None, "You have to define a topology file."
    trj = args["traj"]
    assert trj is not None, "You have to define a trajectory file."
    isomer = args["isomer"]
    assert isomer is not None, "You have to define a isomer type."

    b = int(args["b"])
    e = int(args["e"])
    out += "_{}_{}".format(str(b), str(e))

    if not os.path.exists(out):
        os.mkdir(out)
    else:
        clean_folder(out)

    resid = args["resid"]
    solvent = args["solvent"]
    samples_number = args["NSamples"]

    print("File top:", top)
    print("File trj:", trj)
    print("Time >=", b, "ps")
    print("Time <", e, "ps")
    print("Resid:", resid)

    # Reads the system trajectory
    t = md.load(trj, top=top)

    traj = t[b:e]
    frames_sample = np.random.choice(range(len(traj)), samples_number, replace=False)

    # Extracts the system topology
    top = traj.topology
    table, bonds = top.to_dataframe()

    # The indices of the atoms of interest are selected.
    iat_res = top.select("resid %s" % resid)

    # for frame to gaussian input
    for i in frames_sample:
        atoms = table.loc[iat_res, "element"].values
        coord = traj.xyz[:, iat_res][i] * 10.0  # to angstroms
    
        mol = Molecule(atoms, coord)
    
        # save file gaussian
        GaussianInput(
            mol,
            charge=0,
            spin_multiplicity=1,
            title="Azobenzene %s Sampling in %s - sample %d" % (isomer, solvent, i),
            functional="Cam-B3LYP",
            basis_set="6-311+g(d,p)",
            route_parameters={
                "TD": "(NStates=6)",
                "SCRF": "(Solvent=TetraHydroFuran)"
            },
            link0_parameters={
                "%mem": "8GB",
                "%nprocshared": "12"
            }
        ).write_file("%s/azo%s_%s_%005d.com" % (out, isomer[0].upper(), solvent, i))


if __name__ == '__main__':
    main()
