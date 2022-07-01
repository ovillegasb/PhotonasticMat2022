"""Module dedicated to the analysis of trajectories generated by STAMP."""

import re
import time
import numpy as np
import pandas as pd
from scipy.constants import N_A

""" Regular expression that extracts matrix XYZ """
atoms = re.compile(r"""
        ^\s+
        (?P<atsb>[A-Za-z]+\d?\d?)\s+      # Atom name.
        (?P<x>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for X.
        (?P<y>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Y.
        (?P<z>[+-]?\d+\.\d+)\s+           # Orthogonal coordinates for Z.
        """, re.X)


def read_fatomes(file):
    """Read Fatomes file."""
    natypes = 0
    atomsM = {}
    xyz = []
    connects = dict()
    with open(file, "r") as FATM:
        for line in FATM:
            if "*" == line[0]:
                # ignore lines with the * symbol
                continue

            elif "NbTypesAtomes" in line:
                line = line.split()
                natypes += int(line[1])
                continue

            elif "nom" in line:
                line = line.split()
                nom = line[1]
                continue

            elif "masse" in line:
                line = line.split()
                mass = np.float64(line[1])
                atomsM[nom] = mass
                continue

            elif "maille_long" in line:
                line = line.split()
                box = np.array(line[1:4]).astype(np.float64)
                continue

            elif "PositionDesAtomesCart" in line:
                Natoms = int(FATM.readline())
                continue

            elif atoms.match(line):
                m = atoms.match(line)
                xyz.append(m.groupdict())

            elif "Zmatrice" in line:
                N = int(FATM.readline())
                print("N conectivity:", N)
                for _ in range(N):
                    zline = FATM.readline()
                    zline = zline.split()
                    zline = [int(i) for i in zline]
                    connects[zline[0]] = zline[1:]

    print("Number of atoms in XYZ matrix:", Natoms)

    print("ATOMS types and Mass [Kg/mol]")
    print(atomsM)

    print("Box dimensions [angs]:")
    print(box)

    tabXYZ = pd.DataFrame(xyz)

    tabXYZ = tabXYZ.astype({
        "x": np.float64,
        "y": np.float64,
        "z": np.float64
    })

    tabXYZ["mass"] = tabXYZ["atsb"].apply(lambda x: atomsM[x])

    return tabXYZ, box, connects


def traj_analysis(ndx_mol, top, traj, box, connectivity, GyrationTensor):
    """
    Analyze properties during a simulation.

    Parameters
    ----------
    ndx_mol : dict
        Dictionary with the indexes of each molecule to be analyzed.

    top : DataFrame
        File with the system topology.

    traj : list(DataFrame)
        Defines the trajectory of the system in a list of Dataframes.

    box : numpy.array (1x3)
        Vector box.

    connectivity : molcraft.structure.connectivity
        System connectivity.

    GyrationTensor : class
        Class that builds the information.


    """
    t0 = time.time()
    # Processing bar
    Ntotal = len(traj)
    dframe = int(round(Ntotal / 100))
    state = dict()
    for p in range(0, 100):
        bar = int(p * 20 / 100)
        state[p * dframe] = [
            f"{p*dframe*100/Ntotal:.2f} %",
            bar*"*"+(20-bar)*"_"]

    print("Starting analysis:")
    data = []
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

            G = GyrationTensor(coord, masses, box, pbc=False)
            data.append({
                "frame": i,
                "idx": mol,
                "Nres": int(ndx_mol[mol]["Natoms"] / 10),
                "Rg": G.iso_w_rg,
                "k2": G.shape_anisotropy,
                "dmax": G.max_distance
            })

        if i in state:
            print(f"{state[i][0]} | {state[i][1]} | {time.time()-t0}")

    tf = time.time()
    print(f"Total time: {tf-t0:.2f} sec")
    return pd.DataFrame(data)


def load_data(file, t="NVT"):
    """
    Load data from Stamp.dat.

    Parameters:
    -----------
    file : str
        Output file from STAMP.

    t : str
        Simulation type: NVT or NPT.

    """
    names = {
        "NVT": [
            "I",
            "Etot", "Epot", "Epot_intra", "Epot_inter", "Ekin",
            "T", "Tx", "Ty", "Tz",
            "P", "Px", "Py", "Pz",
            "Vx", "Vy", "Vz",
            "D", "cpu"],
        "NPT": [
            "I",
            "Etot", "Epot", "Epot_intra", "Epot_inter", "Ekin",
            "T", "Tx", "Ty", "Tz",
            "P", "Px", "Py", "Pz",
            "Cmvx", "Cmvy", "Cmvz",
            "Lx", "Ly", "Lz",
            "Vx", "Vy", "Vz",
            "Dens",
            "D", "cpu"
        ]}
    data = pd.read_csv(
        file,
        sep=r"\s+",
        header=None,
        names=names[t],
        comment="#"
                      )
    data["Etot"] = data["Etot"] * N_A / 1000  # to kJ/mol
    return data