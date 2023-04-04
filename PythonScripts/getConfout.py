"""Generates a final gromacs configuration file."""

from stamptools.stamp import STAMP
import numpy as np
from multiprocessing import Pool
from stamptools.analysis import center_of_mass, minImagenC
import sys


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


def get_nopbc_mol(mol_ndx, confout, connectivity, box):
    """Return the structure complet fo a molecule."""
    mol_xyz = confout.loc[mol_ndx, :]
    mol_conn = connectivity.sub_connect(mol_ndx)
    mol_conn.update_coordinates(mol_xyz)
    mol_conn.noPBC(box, center=np.zeros(3))
    mol_conn = mol_conn.reset_nodes()
    mol_conn.simple_at_symbols(add_mass=True)
    return mol_conn.get_df(), mol_ndx


print("Specify the donnees file.")
try:
    donnees = sys.argv[1]
except IndexError:
    print("A Donnees.in file must be specified.")
    print("\tgetConfout.py file.in")
    exit()

# Load System
system = STAMP(donnees=donnees)

resid_PC = 0
atoms_per_mol = system.atoms_per_mol
for resid in atoms_per_mol:
    if resid == resid_PC:
        atoms_per_mol[resid]["resname"] = "MOL"
    elif len(atoms_per_mol[resid]["index"]) == 1:
        atoms_per_mol[resid]["resname"] = "MOL"
    else:
        atoms_per_mol[resid]["resname"] = "MOL"

RESinfo = {}
for resid in atoms_per_mol:
    resname = atoms_per_mol[resid]["resname"]
    index = atoms_per_mol[resid]["index"]
    for i in index:
        RESinfo[i] = {"resid": resid, "resname": resname}

# get trajectory
confout = system.get_traj(b=-1)[0]

# BOXs
box = system.box_in_frame[0][0:3]

print("Saving GRO files")
print("Generating a gro file without PBC")

confout["resid"] = confout.index.map(lambda x: RESinfo[x]["resid"])
confout["resname"] = confout.index.map(lambda x: RESinfo[x]["resname"])

connectivity = system.connectivity

pc_xyz = confout[confout["resid"] == 0]
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

confout["resid"] += 1

save_gro(confout, "confout", box)
print("\nSaved gro file: \033[1;36m%s\033[m writed\n" % "confout.gro")