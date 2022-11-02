#!/bin/env python


"""Script que extrae 100 estructuras."""

import os
import glob
from numpy import random
import numpy as np
import pandas as pd
from stamptools.analysis import read_fatomes, PBC_distance, center_of_mass, translate_to
from molcraft import structure


def Collect_configurations(isomer, replicas, path):
    N_total = 100
    N_choice = int(round(N_total / len(replicas)))
    frames_sample = {}
    print("N frame per replica:", N_choice)
    for i in replicas:
        path_rep = path.format(isomer, i)
        print(path_rep)
        frames_sample[path_rep] = []
        if os.path.exists(path_rep):
            frames = glob.glob("{}/XYZ/*.xyz".format(path_rep))
            N_frames = len(frames)
            # 75 % after
            frame_i = int(round(N_frames * 0.75))
            select_frames = frames[frame_i:N_frames]
            print("N frames select:", len(select_frames))
            # print(random.choices(select_frames, k=N_choice))
            
            frames_sample[path_rep] += list(random.choice(select_frames, size=N_choice, replace=False))
        else:
            print("Folder: {} dont exist.".format(path_rep))
            break

    ntot = 0
    for k in frames_sample:
        ntot += len(frames_sample[k])

    print("N total frames sampling:", ntot)
                  
    return frames_sample
    

def mol_traj_cut_distance(
    connectivity,
    traj,
    atoms_per_mol,
    box,
    top,
    ref=0,
    rcutoff=0.3,
    out_folder="mol_distances",
    comment=""
):
    """Extract the structure of mol around a reference mol using atom-atom distance."""
    print(f"Extract the structure of molecules around a reference mol, resid {ref},")
    print(f"using atom-atom from a distance cut: {rcutoff:.2f} nm.")

    out_folder = out_folder + f"_{rcutoff}"
    rcutoff *= 10.
    masses = top.loc[atoms_per_mol[ref]["index"], "mass"].values
    
    # Create folder where the files will be saveds
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    for frame in traj:
        coord = traj[frame]
        name = "mol%s_%04d" % (comment, frame)
        atoms_ref = atoms_per_mol[ref]["index"]
        mol_conn = connectivity.sub_connect(atoms_ref)
        mol_ref = coord.loc[atoms_ref, :]
        mol_xyz = mol_ref.loc[:, ["x", "y", "z"]].values

        mols_around = []
        for j in atoms_per_mol:
            if j != ref:
                # print(j)
                atoms_near = atoms_per_mol[j]["index"]
                mol_conn_near = connectivity.sub_connect(atoms_near)
                mol_ref_near = coord.loc[atoms_near, :]
                mol_xyz_near = mol_ref_near.loc[:, ["x", "y", "z"]].values
                # print(mol_xyz_near)
                ##### DISTANCE Search
                for v_ref in mol_xyz:
                    new_mol = False
                    for v_near in mol_xyz_near:
                        r = PBC_distance(v_ref, v_near, box)
                        if r <= rcutoff:
                            mols_around.append(j)
                            new_mol = True
                            break
                    if new_mol:
                        break
                #####

        # update coordinates
        mol_conn.update_coordinates(mol_ref)
        # remove PBC
        mol_conn.noPBC(box, center=np.zeros(3))
        mol_coord = mol_conn.get_df()
        co = mol_coord.loc[:, ["x", "y", "z"]].values
        center = center_of_mass(co, masses)
        resids = [ref] + mols_around
        ncoords = []
        for res in resids:
            atoms = atoms_per_mol[res]["index"]
            mol_conn = connectivity.sub_connect(atoms)
            df = coord.loc[atoms, :]
            # Translate to a center of reference
            df.loc[:, ["x", "y", "z"]] = translate_to(df.loc[:, ["x", "y", "z"]].values, center, box)
            # update coordinates
            mol_conn.update_coordinates(df)
            # remove PBC
            mol_conn.noPBC(box, center=np.zeros(3))
            # Reset index and symbols, and add mass
            mol_conn = mol_conn.reset_nodes()
            mol_conn.simple_at_symbols()

            # Search and add hydrogen to vacant atoms
            mol_conn.add_hydrogen(box, type_add="terminal")
            new_mol_xyz = mol_conn.get_df()
            ncoords.append(new_mol_xyz)
            
        ncoords = pd.concat(ncoords, ignore_index=True)
        structure.save_xyz(ncoords, name=f"{out_folder}/{name}")


def gen_files_xyz(frames, home, out="xyz_center"):
    sim = 0
    nframes_analysis = 0
    for r in frames:
        print(r)
    
        # Reading FAtome
        topology, box, connects = read_fatomes(f"{r}/FAtomes.inp_6_prod")
        
        # Connectivity
        connectivity = structure.connectivity()
        connectivity.define_atoms(topology)
        connectivity.read_dict(connects)
        
        atoms_per_mol = connectivity.atomsMOL
        
        # trajectory
        XYZs = frames[r]
        nframes_analysis += len(XYZs)
        traj = {}
        for file in XYZs:
            print(file)
            n_frame = file.split("/")[-1].split("_")[-1].split(".")[0]
            print(n_frame)
            traj[int(n_frame)] = structure.load_xyz(file, warning=False)

        # Gen configuration rcutoff 0.3 nm
        # print(traj)
        
        mol_traj_cut_distance(
            connectivity,
            traj,
            atoms_per_mol,
            box,
            topology,
            out_folder=f"{home}/{out}",
            comment=f"_{sim}"
        )
        sim += 1

    print("N total frames analized:", nframes_analysis)


def main():
    folder = "../azo{}_procedure/6_prod_{}"
    outfolder = "../ONIOM_gaus/azoOC_sampling/sampling_100_method3"
    frames = Collect_configurations("C", ["0_long", 1, 2, 3, 4], folder)

    gen_files_xyz(frames, outfolder)


if __name__ == '__main__':
    main()
