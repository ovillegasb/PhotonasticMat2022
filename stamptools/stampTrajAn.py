"""Program to automate the analysis of trajectories in STAMP."""

import sys
import argparse
import glob
import time

# Adding some modules path - modules extras
sys.path.append("/home/ovillegas/GITPROYECTS/molcraft/")
sys.path.append("/home/ovillegas/GITPROYECTS/PhotonasticMat/")

from molcraft import structure  # noqa: E402
from molcraft import clusters  # noqa: E402
from analysis import read_fatomes  # noqa: E402


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="StampTrajAn",
        usage="%(prog)s FAtomes.in [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    # FAtomes.in file
    parser.add_argument(
        "fatoms",
        help="FAtomes.in file",
        default="",
        type=str
    )

    return vars(parser.parse_args())


# start time
t0 = time.time()
args = options()
fatoms = args["fatoms"]
print(f"FAtomes file: {fatoms}")

# Reading FAtome
top, box, connects = read_fatomes(fatoms)

# You get a list of all xyz files to be analyzed.
# XYZs = glob.glob("XYZ/PasDeCalcul_Iteration_00200000.xyz")
XYZs = sorted(glob.glob("XYZ/*.xyz"))

# Builds the system trajectory
traj = list()

for file in XYZs:
    xyz = structure.load_xyz(file, warning=False)
    traj.append(xyz)

# Number of frames read
Nframes = len(traj)
print(f"Number of frames: {Nframes}")

# Brings system connectivity
connectivity = structure.connectivity()
connectivity.define_atoms(top)
connectivity.read_dict(connects)

# atoms per polymer
ndx_mol = connectivity.atomsMOL

# Processing bar
dframe = int(round(Nframes / 100))
state = dict()
for p in range(0, 100):
    bar = int(p * 20 / 100)
    progress = bar*"=" + (20-bar)*" "
    state[p * dframe] = f"{p*dframe*100/Nframes:5.2f} % |{progress}|"
t1 = time.time()
print(f"time {t1-t0:.2f} s")

out = open("polymers.csv", "w")
out.write(",frame,idx,Nres,Rg,k2,dmax\n")
print("\nStarting analysis:")
index = 0
for i, frame in enumerate(traj):
    for mol in ndx_mol:
        masses = top.loc[ndx_mol[mol]["index"], "mass"].values
        dfcoord = frame.loc[ndx_mol[mol]["index"], :]

        # Connectivity in the molecule
        connect = connectivity.sub_connect(ndx_mol[mol]["index"])

        # update coordinates
        connect.update_coordinates(dfcoord)

        # remove PBC
        connect.noPBC(box)

        newdfcoord = connect.get_df()
        coord = newdfcoord.loc[:, ["x", "y", "z"]].values

        G = clusters.GyrationTensor(coord, masses, box, pbc=False)
        line = ""
        line += f"{index},"
        line += f"{i},"
        line += f"{mol},"
        line += "{},".format(int(ndx_mol[mol]["Natoms"] / 10))
        line += f"{G.iso_w_rg:.2f},"
        line += f"{G.shape_anisotropy:.3f},"
        line += f"{G.max_distance:.2f},"
        line += "\n"

        out.write(line)
        index += 1

    try:
        print(f"{state[i]} - {time.time()-t1:.2f} s")
    except KeyError:
        pass

out.close()

tf = time.time()
print(f"Analysis time: {tf-t1:.2f} s")
