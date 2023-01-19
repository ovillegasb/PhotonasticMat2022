"""Submodule dedicated to the analysis of trajectories generated by GROMACS."""

import time
import mdtraj as md


def read_xtc(**kwargs):
    """Read the trajectrory using mdtraj module from the xtc file."""
    # set needed variables:
    trajfile = kwargs["xtc"]
    top = kwargs["top"]
    interval = kwargs["interval"]
    begin = kwargs["b"]
    end = kwargs["e"]

    print("Read trajectory file: ", trajfile)
    print("Read topology from:   ", top)

    t1 = time.time()
    print("Start reading at", time.strftime("%H:%M:%S", time.localtime(t1)))
    print("  first frame: ", begin)
    print("  last frame:  ", end if end > 0 else "last")
    print("  step:        ", interval, "\n")
    
    trajectory = md.load(trajfile, top=top)
    t2 = time.time() - t1
    print("Done in %.0f s" % t2)

    print("  # atoms:           ", trajectory.n_atoms)
    print("  # frames total:    ", trajectory.n_frames)
    if trajectory.n_frames > 1:
        trajectory = trajectory[begin: end: interval]
    print("  # frames selected: ", trajectory.n_frames, "\n")

    return trajectory
