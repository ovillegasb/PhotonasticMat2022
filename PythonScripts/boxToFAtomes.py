#!/usr/bin/env python
# -*- coding=utf-8 -*-

import argparse
import mdtraj as md
from molcraft.structure import MOL, BULK, Elements
from stamptools.fatomes import FATOME, TOPOL


def options():
    """Generate command line interface."""
    parser = argparse.ArgumentParser(
        prog="boxToFAtomes",
        usage="%(prog)s [-options]",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Enjoy the program!"
    )

    fileinput = parser.add_argument_group(
        """\033[1;36mInitial settings\033[m

Example:

python ~/.gitproyects/PhotonasticMat/PythonScripts/boxToFAtomes.py -box box.gro -s ../../0_files/{azoOT.pdb,btdn.pdb} -dummyPC

"""
    )

    fileinput.add_argument(
        "-box", "--box",
        help="System box in gro format.",
        type=str
    )

    fileinput.add_argument(
        "-f", "--fatomes",
        help="The system is read from a FAtomes.",
        type=str,
        default=None
    )

    fileinput.add_argument(
        "-s", "--topols",
        help="Input structure, system topology: pdb format. It must be defined\
in the same order as the file is written.",
        type=str,
        nargs="+",
        default=[]
    )

    fileinput.add_argument(
        "-dummyPC", "--dummyPC",
        help="Generates a FAtome with DU-type photochrome atoms.",
        action="store_true"
    )

    fileinput.add_argument(
        "--addPC",
        help="Adds the FF parameters of the photochrome to the fatomes.",
        action="store_true"
    )

    fileinput.add_argument(
        "--add_OH",
        help="Reads a fatome and adds an OH group to vacant terminal carbons.",
        action="store_true"
    )

    return vars(parser.parse_args())


def main():
    """Run program."""
    args = options()
    print(args)

    # extract parameters
    box = args["box"]
    fatomes = args["fatomes"]

    # search connectivities files
    top_files = args["topols"]
    print("Topology files:", top_files)

    if box is None:
        print("No box has been taken")
        print("It will be assumed that it will be read from a fatomes.")
        assert fatomes is not None, "HEY! you haven't taken any FAtomes."

        fatomes = TOPOL(fatomes)
        print("FAtomes file:", fatomes.file)

        if args["addPC"]:
            assert len(top_files) > 0, "To add the PC you need a topology file."
            fatomes.activate_PC(PCtopol=top_files[0])

        if args["add_OH"]:
            fatomes.complete_with_OH()

        fatomes.export_xyz("systems.xyz")
        # File
        output = args["fatomes"].split("/")[-1].replace(".in", "_modif.in")
        fatomes.save_fatomes(
            name=output
        )

    else:
        t = md.load(box, top=box)

        # Extracts the system topology
        top = t.topology
        table, bonds = top.to_dataframe()
        dfatoms = table.copy()
        dfatoms["atsb"] = dfatoms["element"]
        dfatoms["x"] = t.xyz[0][:, 0]
        dfatoms["y"] = t.xyz[0][:, 1]
        dfatoms["z"] = t.xyz[0][:, 2]
        dfatoms["mass"] = dfatoms["atsb"].apply(lambda at: Elements[at]["mass"])
        dfatoms["charge"] = 0.0
        print(dfatoms)

        resnames = list(table["resName"].unique())
        print("ResName:", " ".join(resnames))

        nres = len(resnames)
        print("Number of restype:", nres)
        if args["dummyPC"] and nres > 1:
            coordPC = dfatoms[dfatoms["resName"] == resnames[0]].copy()
            coordPC["name"] = "DU"
            dfatoms[dfatoms["resName"] == resnames[0]] = coordPC

        if not len(top_files) == nres:
            print("You must use the same number of topology files as residue files.")
            exit()

        # Se lee los archivos con la conectividad
        SystemConnectivity = {}
        for i, res in enumerate(resnames):
            # print("RESNAME:", res, i)
            mol = MOL(file=top_files[i], res=res, connectivity=True)
            connect = mol.connect.atoms_map
            # print(connect)
            # print(mol)
            rescoords = table[table["resName"] == res]
            # print(rescoords)
            res_j = list(rescoords["resSeq"].unique())
            # print("res i:", res_j)
            for n in res_j:
                # print("RES:", n)
                coord = rescoords[rescoords["resSeq"] == n]
                factor = 0
                for j, at_j in enumerate(coord.index):
                    # index = i + n - 1
                    if j == 0:
                        factor = at_j
                
                    # print(j, connect[j], "--->", at_j, [val + factor for val in connect[j]])
                    SystemConnectivity[at_j] = [val + factor for val in connect[j]]

        # Creates the object BULK
        bulk = BULK(box=t.unitcell_lengths[0])
        bulk.read_connectivity(dfatoms, SystemConnectivity)

        # Creates the FAtomes object and reads the created BULK object.
        fatomes = FATOME()
        fatomes.read_connect(bulk)

        # Treatement
        fatomes.get_atomstypes()

        # Save file
        fatomes.save_file()


if __name__ == '__main__':
    main()
