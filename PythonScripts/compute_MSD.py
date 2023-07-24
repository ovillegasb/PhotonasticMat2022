#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""Compute the MSD using the MDAnlysis module."""

import argparse
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import pandas as pd


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="compute_MSD",
        usage="%(prog)s [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    fileinput = parser.add_argument_group(
        "\033[1;36mInitial settings\033[m"
    )

    fileinput.add_argument(
        "-f", "--traj",
        help="Input trajectory or single configuration: xtc, gro, xyz.",
        type=str,
        default=None
    )

    fileinput.add_argument(
        "-s", "--topol",
        help="Input structure, system topology: gro.",
        type=str,
        default=None
    )

    fileinput.add_argument(
        "-dt", "--dt",
        help="Time step in ps",
        type=float,
        default=1.0
    )

    fileinput.add_argument(
        "-o", "--out",
        help="Output file.",
        type=str,
        default="msd.csv"
    )

    fileinput.add_argument(
        "--select",
        help="Selection type VMD",
        nargs="+",
        default=["index", "0" "to" "25"]
    )

    return vars(parser.parse_args())


def main():
    """Run program."""
    args = options()

    out = args["out"]
    dt = args["dt"]

    top = args["topol"]
    assert top is not None, "You have to define a topology file."
    traj = args["traj"]
    assert traj is not None, "You have to define a trajectory file."

    universe = mda.Universe(top, traj, dt=dt)

    # Select Photochrome
    select = " ".join(args["select"])  # "index 0 to 25"

    # Class MSD
    MSD = msd.EinsteinMSD(universe, select=select, msd_type='xyz', fft=True)

    # compute MSD
    MSD.run()

    # get msd
    msd_data = MSD.results.timeseries
    
    nframes = MSD.n_frames
    lagtimes = np.arange(nframes)*dt  # make the lag-time axis

    dfmsd = pd.DataFrame({"time": lagtimes, "msd": msd_data})
    dfmsd.to_csv(out, float_format="%e", index=False)


if __name__ == '__main__':
    main()
