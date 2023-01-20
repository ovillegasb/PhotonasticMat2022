"""RDF analysis all atoms between systems."""

import glob
import numpy as np
from stamptools.analysis import read_fatomes, PBC_distance
from molcraft.structure import connectivity, load_xyz
import sys
from scipy.spatial.distance import cdist
import time
import pandas as pd

# Folder and systems
systems = ["0", "1", "2", "3", "4", "0_long"]
isomer = "azoC_procedure"
home = "/home/ovillegas/.bettyboop"
fatomes = f"{home}/{isomer}/6_prod_0/FAtomes.inp_6_prod"
name = "RDF_all_test"

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
    if i == 50:
        break

end = time.perf_counter()

print("frame_distances computed: " + str(round(end - start, 2)) + "s")
print("Nomber of frames added:", len(frame_distances))
print(
    "Memory in frame_distances:",
    sys.getsizeof(frame_distances),
    "bytes"
)

# RDF
rdf_start = time.perf_counter()

n_centers = len(atoms_ref)
n_centers_near = len(atoms_env)

vol_per_sphere = vol / n_centers_near
vshell = 4 * np.pi * ((binwidth + bins)**3 - bins**3) / 3


g_r = np.zeros(len(bins))
for frame in frame_distances:
    for i, atom in enumerate(frame_distances[frame]):
        indexs = np.int64(atom * 0.1 / binwidth)
        indexs = indexs[indexs < len(bins)]
        for n in indexs:
            g_r[n] += 1

n_frames = len(frame_distances.keys())
g_r_norm = g_r * vol_per_sphere / vshell / n_frames / n_centers

rdf_end = time.perf_counter()
print("RDF computed: " + str(round(rdf_end - rdf_start, 2)) + "s")


# return g_r_norm, bins
RDF = pd.DataFrame({
    "g_r": g_r_norm,
    "r": bins
})

file = f"{name}.csv"
RDF.to_csv(file, float_format="%.6f", index=None)
print(f"file {file} saved.", end=" - ")
print("")
