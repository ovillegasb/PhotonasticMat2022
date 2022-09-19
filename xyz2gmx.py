#!/bin/env python

"""Program created to convert .xyz files to .gro files."""

import argparse
import json
import glob
import numpy as np
from molcraft import structure
from stamptools.analysis import read_fatomes


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="xys2gro",
        usage="%(prog)s [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    fileinput = parser.add_argument_group(
        "\033[1;36mInput files\033[m")

    # FAtomes.in file
    fileinput.add_argument(
        "-f", "--fatomes",
        help="FAtomes.in file",
        default="",
        type=str
    )

    # XYZ files
    fileinput.add_argument(
        "-xyz", "--xyz",
        help="XYZ files",
        nargs="+",
        default=None
        )


def load_conectivity(file):
    """Open JSON file."""
    with open(file) as json_file:
        data = json.load(json_file)

    return data


def save_gro(name, description, coord, box):
    """Write coordinates to file gro."""
    nat = len(coord)

    GRO = open('%s.gro' % name, 'w', encoding='utf-8')
    GRO.write('%s\n' % description)
    GRO.write('%5d\n' % nat)
    for i in coord.index:
        GRO.write('{:>8}{:>7}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(
            coord.loc[i, 'resid'] + coord.loc[i, 'resname'],
            coord.loc[i, 'atsb'].upper(),
            i + 1,
            coord.loc[i, 'x'] * 0.1,
            coord.loc[i, 'y'] * 0.1,
            coord.loc[i, 'z'] * 0.1)
        )
    GRO.write('   {:.5f}   {:.5f}   {:.5f}\n'.format(
        box[0] * 0.1,
        box[1] * 0.1,
        box[2] * 0.1))

    GRO.close()
    print("file %s.gro writed" % name)


def main():
    """Program in terminal."""
    args = options()
    print(args)
    exit()
    # Reading FAtome
    _, box, conn = read_fatomes("FAtomes.inp_6_prod")

    print(conn)

    # Load conectivity
    connect = load_conectivity("polysystem.json")

    # Load time per frame
    t_frame = np.loadtxt("t_frame.dat")

    # Load XYZ trajectory
    XYZ = glob.glob("XYZ/PasDeCalcul_*")
    traj = list()

    for file in XYZ:
        xyz = structure.load_xyz(file, warning=False)
        traj.append(xyz)

    # Gen resid and resname list
    resid = list()
    resname = list()
    for rid in connect:
        idx = connect[rid]["index"].split()
        res = "POL"
        if rid == "0":
            res = "AZO"
        for _ in idx:
            resid.append(rid)
            resname.append(res)

    for i, coord in enumerate(traj):
        coord["resname"] = resname
        coord["resid"] = resid
        coord["x"] += box[0]/2
        coord["y"] += box[1]/2
        coord["z"] += box[2]/2
        description = "MD PC - Polymer, t=%.3f" % t_frame[i]
        print(description)
        save_gro("GRO/frame_%005d" % i, description, coord, box)
        pass


if __name__ == '__main__':
    main()
