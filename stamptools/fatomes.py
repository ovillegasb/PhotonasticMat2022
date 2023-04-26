"""Module to define the FATOME class."""

from stamptools.analysis import decoTime, change_atsb
from stamptools.analysis import center_of_mass, translate_to
from molcraft import structure
import itertools as it
import numpy as np
import pandas as pd
import networkx as nx
from datetime import datetime
import re
import os


# dat location
location = os.path.dirname(os.path.realpath(__file__))
cell_unit = os.path.join(location, "oplsaa.dat")


# Some definitions
POTENTIAL_par_OPLS = {
    "OH": "LJ sigma 3.12 ang epsilon 0.066 kcal/mol rc 12.5 ang",
    "HO": "LJ sigma 0.00 ang epsilon 0.000 kcal/mol rc 12.5 ang",
    "CT": "LJ sigma 3.50 ang epsilon 0.066 kcal/mol rc 12.5 ang",
    "HT": "LJ sigma 2.50 ang epsilon 0.030 kcal/mol rc 12.5 ang"
}

BONDS_par_OPLS = {
    ("HO", "OH"): "bond_opls     HO  OH     0.945 ang   553.0 kcal/mol/ang2   0.0 -  0.0 -",
    ("CT", "OH"): "bond_opls     CT  OH     1.410 ang   320.0 kcal/mol/ang2   0.0 -  0.0 -"
}

ANGLES_par_OPLS = {
    ("CT", "OH", "OH"): "angle_opls     CT  OH  HO     108.5 degre  55.0  kcal/mol   0.0 -  0.0 -",
    ("CM", "CT", "OH"): "angle_opls     CM  CT  OH     109.5 degre  50.0  kcal/mol   0.0 -  0.0 -",
    ("HT", "CT", "OH"): "angle_opls     HT  CT  OH     109.5 degre  35.0  kcal/mol   0.0 -  0.0 -"
}

# https://pubs-acs-org.inc.bib.cnrs.fr/doi/suppl/10.1021/ja9621760/suppl_file/ja11225.pdf
DIHEDRALS_par_OPLS = {
    ("CM", "CT", "OH", "HO"): "torsion_opls     CM  CT  OH  HO  -0.356 kcal/mol -0.174  kcal/mol   0.492  kcal/mol",  # amberTools 051
    ("HT", "CT", "OH", "HO"): "torsion_opls     HT  CT  OH  HO   0.000 kcal/mol  0.000  kcal/mol   0.450  kcal/mol",  # oplsaa.dat
    ("CM", "CM", "CT", "OH"): "torsion_opls     CM  CM  CT  OH   1.711 kcal/mol -0.500  kcal/mol   0.663  kcal/mol",
    ("HM", "CM", "CT", "OH"): "torsion_opls     HM  CM  CT  OH   0.000 kcal/mol  0.000  kcal/mol   0.468  kcal/mol"
}

Keywords = {}


""" Regular expression for FAtomes """
noms = re.compile(r"""
    nom\s+(?P<nom>\w+)
    """, re.X)

masses = re.compile(r"""
    masse\s+(?P<mass>\d+\.\w+[+-]?\d+)
    """, re.X)

charges = re.compile(r"""
    charge\s+(?P<charge>[+-]?\d+\.\d+\w?[+-]?\d*)
    """, re.X)

""" Regular expression that extracts matrix XYZ """
atoms = re.compile(r"""
        ^\s*
        (?P<atypes>[A-Za-z]+\d?\d?)\s+            # Atom name.
        (?P<x>[+-]?\d+\.\d+\w?[+-]?\d*)\s+      # Orthogonal coordinates for X.
        (?P<y>[+-]?\d+\.\d+\w?[+-]?\d*)\s+      # Orthogonal coordinates for Y.
        (?P<z>[+-]?\d+\.\d+\w?[+-]?\d*)\s+      # Orthogonal coordinates for Z.
        """, re.X)


class TOPOL:
    """Object describing the system topology."""

    def __init__(self, fatomes):
        """Initialize the class by reading the indicated Fatomes file."""
        self.file = fatomes
        self._read_fatomes()

        # searching conectivity
        self._get_connectivity()

    @decoTime
    def _read_fatomes(self):
        """Read Fatomes file."""
        print("Reading fatomes file", end=" - ")

        natypes = 0
        atomsM = {}
        xyz = []
        lnom = []
        lmass = []
        lcharges = []
        lcharge_mod = []
        connects = dict()

        # fatome_lines = []

        # ATOMS types
        atoms_types = {}
        infoType = {
            "nom": "",
            "nomXYZ": "",
            "nomFF": "",
            "type": "",
            "masse": "",
            "charge": "",
            "gele": "",
            "structure": "",
            "maille_long": "",
            "maille_angle": "",
            "maille_orient": "",
            "maille_ref": ""
        }
        ntype = 0

        # Intermolecular parameters
        Intermol_potentials = []
        comb_rules = "" 

        # Force field parameters
        ffparms = []

        # Charges modified
        ChargesMOD = []

        # Contributions intramolecular to dispersion
        ContribDispRepIntra = []

        with open(self.file, "r") as FATM:
            nom = ""
            for line in FATM:
                # fatome_lines.append(line)
                if line.startswith("*"):
                    # ignore lines with the * symbol
                    continue

                if noms.match(line):
                    nom = ""
                    m = noms.match(line)
                    lnom.append(m.groupdict())
                    nom += m.groupdict()["nom"]
                    ntype += 1

                if ntype > 0 and ntype <= natypes:
                    break_line = line.split()
                    keywords = break_line[0]
                    if keywords in infoType:
                        atoms_types[ntype - 1][keywords] = " ".join(break_line[1:])

                if masses.match(line):
                    m = masses.match(line)
                    lmass.append(m.groupdict())

                if charges.match(line):
                    m = charges.match(line)
                    lcharges.append((nom, float(m.groupdict()["charge"])))
                    nom = ""

                if "NbTypesAtomes" in line:
                    line = line.split()
                    natypes += int(line[1])
                    atoms_types = {idx: dict(infoType) for idx in range(0, natypes)}
                    continue

                if "maille_long" in line:
                    line = line.split()
                    box = np.array(line[1:4]).astype(np.float64)
                    continue

                if atoms.match(line):
                    m = atoms.match(line)
                    xyz.append(m.groupdict())

                if "ChampDeForces" in line:
                    N = int(FATM.readline())
                    for _ in range(N):
                        lne = FATM.readline()
                        while lne.startswith("*"):
                            lne = FATM.readline()

                        ffparms.append(lne.replace("\n", ""))

                if "Zmatrice" in line:
                    N = int(FATM.readline())
                    # print("N conectivity:", N)
                    for _ in range(N):
                        zline = FATM.readline()
                        zline = zline.split()
                        zline = [int(i) for i in zline]
                        connects[zline[0]] = zline[1:]

                if "Regle_melange" in line:
                    comb_rules = line.replace("\n", "").split()[1]

                if "ModificationChargeDesAtomes" in line:
                    N = int(FATM.readline())
                    for _ in range(N):
                        chline = FATM.readline()
                        ChargesMOD.append(chline.replace("\n", ""))
                        chline = chline.split()
                        lcharge_mod.append({"idx": int(chline[0]), "charge": float(chline[1])})

                if "ContribDispRepIntra" in line:
                    N = int(FATM.readline())
                    for _ in range(N):
                        lne = FATM.readline()
                        while lne.startswith("*"):
                            lne = FATM.readline()

                        ContribDispRepIntra.append(lne.replace("\n", ""))

                if line.startswith("Potentiel"):
                    Intermol_potentials.append(line.replace("Potentiel  ", "").replace("\n", ""))

        lcharges = {key: charge for (key, charge) in lcharges}
        
        def from_lcharges(x):
            """Return partial charge for the atomic index."""
            try:
                return lcharges[x]
            except KeyError:
                return 0.0
        
        lcharge_mod = pd.DataFrame(lcharge_mod)
        atomsM = dict()

        if len(lmass) == len(lnom):
            for i in range(len(lnom)):
                atomsM[lnom[i]["nom"]] = np.float64(lmass[i]["mass"])
                if lnom[i]["nom"][0].upper() == "H":
                    atomsM[lnom[i]["nom"]] = np.float64(1.008e-03)
        else:
            print("\nERROR, no; and masses dont similar")
            exit()

        tabXYZ = pd.DataFrame(xyz)
        tabXYZ["atsb"] = tabXYZ["atypes"].apply(change_atsb)
        tabXYZ = tabXYZ.astype({
            "x": np.float64,
            "y": np.float64,
            "z": np.float64
        })

        try:
            tabXYZ["mass"] = tabXYZ["atypes"].apply(lambda x: atomsM[x])
        except KeyError:
            print("")
            print(atomsM)
            print(tabXYZ.loc[0:26, :])
            print("ERROR")
            exit()

        tabXYZ["charge"] = tabXYZ["atypes"].apply(from_lcharges)

        for i in lcharge_mod.index:
            tabXYZ.loc[lcharge_mod.loc[i, "idx"], "charge"] = lcharge_mod.loc[i, "charge"]

        self.dfatoms = tabXYZ
        self.box = box
        self.connects = connects

        self.atoms_types = atoms_types

        pairsPot = {}
        for pot in Intermol_potentials:
            pot = pot.split()
            iat = int(pot[0])
            jat = int(pot[1])
            itype = atoms_types[iat]["nomFF"]
            jtype = atoms_types[jat]["nomFF"]
            pars = "  ".join(pot[2:])

            if (iat, itype) not in pairsPot:
                pairsPot[(iat, itype)] = {}

            if (jat, jtype) not in pairsPot[(iat, itype)]:
                pairsPot[(iat, itype)][(jat, jtype)] = pars

        self.Intermol_potentials = Intermol_potentials
        self.pairsPot = pairsPot
        self.comb_rules = comb_rules
        self.ffparms = ffparms
        self.ChargesMOD = ChargesMOD
        self.ContribDispRepIntra = ContribDispRepIntra

    @decoTime
    def _get_connectivity(self):
        """Brings system connectivity."""
        print("Searching connectivity", end=" - ")
        conn = structure.connectivity()
        conn.define_atoms(self.dfatoms)
        conn.read_dict(self.connects)
        self.connectivity = conn

    @property
    def atoms_per_mol(self):
        """List atoms per molecule."""
        return self.connectivity.atomsMOL

    def export_xyz(self, file="AtomesCart"):
        """Extract and save the atomic coordinates from the FAtomes."""
        coord = self.dfatoms.copy()
        coord["atsb"] = coord["atsb"].apply(change_atsb)
        structure.save_xyz(coord, file)
        print(f"Saved file: {file}")

    @decoTime
    def complete_with_H(self):
        """Verify and complete molecules missing terminal hydrogens."""
        print("Adding H atom", end=" - ")
        atoms_per_mol = self.atoms_per_mol
        connectivity = self.connectivity

        ncoords = []
        newconnectivity = structure.connectivity()
        conn_list = []
        for resid in atoms_per_mol:
            atoms = atoms_per_mol[resid]["index"]
            # print("resid:", resid)

            # Extracts the connectivity of the residue.
            mol_conn = connectivity.sub_connect(atoms)
            mol_conn.simple_at_symbols()

            # Reset index and symbols, and add mass
            mol_conn = mol_conn.reset_nodes()

            # If the residue has only one atom, it continues to the following.
            if len(atoms) == 1:
                # mol_xyz = mol_conn.get_df()
                # ncoords.append(mol_xyz)
                conn_list.append(mol_conn)
                continue

            # whole molecule
            mol_conn.noPBC(self.box, center=np.zeros(3))

            # Search and add hydrogen to vacant atoms
            mol_conn.add_hydrogen(type_add="terminal-C", mass=1.008e-03, atypes="HT", charge=0.060)

            # Adding PBC
            mol_conn.addPBC(self.box, center=np.zeros(3))
            # mol_xyz = mol_conn.get_df()

            # ncoords.append(mol_xyz)
            conn_list.append(mol_conn)
            
        # Create a new connectivity object with the new information.
        for mol in conn_list:
            newconnectivity.add_residue(mol)

        self.connectivity = newconnectivity
        ncoords = self.connectivity.get_df()

        self.dfatoms = ncoords.copy()

    @decoTime
    def complete_with_OH(self):
        """Verify and complete molecules missing terminal hydrogens."""
        print("Adding OH atom", end=" - ")
        atoms_per_mol = self.atoms_per_mol
        connectivity = self.connectivity

        ncoords = []
        newconnectivity = structure.connectivity()
        conn_list = []

        for resid in atoms_per_mol:
            atoms = atoms_per_mol[resid]["index"]

            # Extracts the connectivity of the residue.
            mol_conn = connectivity.sub_connect(atoms)
            mol_conn.simple_at_symbols()

            # Reset index and symbols, and add mass
            mol_conn = mol_conn.reset_nodes()

            # If the residue has only one atom, it continues to the following.
            if len(atoms) == 1:
                # mol_xyz = mol_conn.get_df()
                # ncoords.append(mol_xyz)
                conn_list.append(mol_conn)
                continue

            # whole molecule
            mol_conn.noPBC(self.box, center=np.zeros(3))

            # Search and add oxygens to vacant carbons
            mol_conn.add_OH(type_add="terminal-C", mass=1.600e-02, atypes="OH", charge=-0.683)
            # Search and add hydrogen to vacant atoms
            mol_conn.add_hydrogen(type_add="OH", mass=1.008e-03, atypes="HO", charge=0.418)
            mol_xyz = mol_conn.get_df()
            structure.save_xyz(mol_xyz, "testOH")
            # Adding PBC
            mol_conn.addPBC(self.box, center=np.zeros(3))

            conn_list.append(mol_conn)
        
        # Create a new connectivity object with the new information.
        for mol in conn_list:
            newconnectivity.add_residue(mol)

        self.connectivity = newconnectivity
        ncoords = self.connectivity.get_df()

        self.dfatoms = ncoords.copy()

    @decoTime
    def center_to(self, ref, ref_type="atoms_ndx"):
        """Centers the system using a reference."""
        print("Centering the system to a reference", end=" - ")
        dfatoms = self.dfatoms.copy()
        if ref_type == "atoms_ndx":
            print("Reference is atoms index", end=" - ")
            dfref = dfatoms.loc[ref, :]
            ref_xyz = dfref.loc[:, ["x", "y", "z"]].values
            center = center_of_mass(ref_xyz, np.ones(ref_xyz.shape[0]))
            dfatoms.loc[:, ["x", "y", "z"]] = translate_to(
                dfatoms.loc[:, ["x", "y", "z"]].values,
                center,
                self.box
            )

            dfatoms["x"] = dfatoms["x"] - center[0] + 0.5 * self.box[0]
            dfatoms["y"] = dfatoms["y"] - center[1] + 0.5 * self.box[1]
            dfatoms["z"] = dfatoms["z"] - center[2] + 0.5 * self.box[2]

            # update coordinates
            self.connectivity.update_coordinates(dfatoms)
            self.connectivity.addPBC(self.box, center=np.zeros(3))
            self.dfatoms = self.connectivity.get_df()

    def _verify_atomTypes(self):
        """Check the atom types with the modified atoms when adding new atoms."""
        new_types = []
        atoms_types = self.atoms_types.copy()
        NatomsTypes = len(atoms_types)
        modified_atoms = self.connectivity.modified_atoms
        if len(modified_atoms) > 0:
            for at in modified_atoms:
                if modified_atoms[at] not in new_types:
                    new_types.append(modified_atoms[at])

            for i, atyp in enumerate(new_types):
                atoms_types[i + NatomsTypes] = atyp

        self.atoms_types = atoms_types
        self.new_types = new_types

    def get_atomstypesPot(self):
        """Check that each type of atom contains an intermolecular van der waals potential."""
        atoms_types = self.atoms_types
        pairsPot = self.pairsPot

        for atom in atoms_types:
            patom = (atom, atoms_types[atom]["nomFF"])
            if patom not in pairsPot:
                iatom = patom[0]
                itype = patom[1]
                if itype in POTENTIAL_par_OPLS:
                    pairsPot[patom] = {}
                    pairsPot[patom][patom] = POTENTIAL_par_OPLS[itype]

                else:
                    print("Parameter not found:", iatom, itype)
                    exit()

        self.pairsPot.update(pairsPot)

    def check_newsFFpar(self):
        """List the force field parameters from the newly added parameters."""
        atypes = []
        pairsPot = self.pairsPot
        for iat, itype in pairsPot:
            atypes.append(itype)

        types_added = []
        for par in self.ffparms:
            if "bond_opls" in par:
                types_added.append(tuple(par.split()[1:3]))

        new_types = [at["nomFF"] for at in self.new_types]
        for i, j in it.combinations(new_types, 2):
            if (i, j) in types_added or (j, i) in types_added:
                pass
            elif i[0].upper() == "H" and j[0].upper() == "H":
                pass
            elif i == "HT" and j == "OH":
                pass
            elif j == "HT" and i == "OH":
                pass
            elif i == "CT" and j == "HO":
                pass
            elif j == "CT" and i == "HO":
                pass
            else:
                try:
                    self.ffparms.append(BONDS_par_OPLS[(i, j)])
                    types_added.append((i, j))
                except KeyError:
                    self.ffparms.append(BONDS_par_OPLS[(j, i)])
                    types_added.append((j, i))

                # angles
                for ang in ANGLES_par_OPLS:
                    if i in ang or j in ang:
                        if ANGLES_par_OPLS[ang] not in self.ffparms:
                            self.ffparms.append(ANGLES_par_OPLS[ang])

                # torsion
                for dih in DIHEDRALS_par_OPLS:
                    if i in dih or j in dih:
                        if DIHEDRALS_par_OPLS[dih] not in self.ffparms:
                            self.ffparms.append(DIHEDRALS_par_OPLS[dih])

    @decoTime
    def save_fatomes(self, name=None):
        """Save the fatomes in a file."""
        print("Writing FAtomes file:", end=" - ")
        now = datetime.now()
        now = now.strftime("%d/%m/%y %H:%M:%S")
        header = f"""*
*#####################################################################*
*                       Fichier issu de Stamp                         *
*                        (Iteration   302500)                         *
*                                                                     *
*        Date: {now}                                      *
*        Created using STAMPtools                                     *
*#####################################################################*
*
* ======================= 
* Bibliotheque des atomes 
* =======================
*
"""

        out = ""
        if name:
            out += name.replace(".in", "") + ".in"
        else:
            out += self.file.replace(".in", "") + "_modif.in"

        lines = ""
        # header
        lines += header

        self._verify_atomTypes()
        atoms_types = self.atoms_types
        # Atoms types
        # ==============================
        lines += "NbTypesAtomes %d\n" % len(atoms_types)

        for i in atoms_types:
            type_lines = f"* Description du type atomique    {i}\n"
            for key in atoms_types[i]:
                if atoms_types[i][key] != "":
                    type_lines += "{} {}\n".format(key, atoms_types[i][key])
            
            lines += type_lines

        # Intermoleculars potentials
        # ==============================
        lines += """*
* ============================ 
* Potentiels intermoleculaires 
* ============================ 
"""
        self.get_atomstypesPot()
        pairsPot = self.pairsPot
        for patom1 in pairsPot:
            for patom2 in pairsPot[patom1]:
                pars = pairsPot[patom1][patom2]
                iat = patom1[0]
                jat = patom2[0]
                lines += f"Potentiel  {iat}  {jat}  {pars}\n"

        if self.comb_rules != "":
            lines += "Regle_melange %s\n" % self.comb_rules

        # Forcefield parameters
        # ==============================
        lines += """*
* =============== 
* Champ de forces 
* =============== 
ChampDeForces
"""     
        self.check_newsFFpar()
        ffparms = self.ffparms
        lines += "%d\n" % len(ffparms)
        for p in ffparms:
            lines += "%s\n" % p

        # Atomics positions
        # ==============================
        lines += """*
* =================== 
* Positions atomiques 
* =================== 
PositionDesAtomesCart ang
"""
        dfatoms = self.dfatoms.copy()
        lines += "%d\n" % len(dfatoms)
        for i in dfatoms.index:
            line_xyz = (dfatoms.atypes[i], dfatoms.x[i], dfatoms.y[i], dfatoms.z[i])
            lines += '%5s%15.6f%15.6f%15.6f\n' % line_xyz

        # System connectivity
        # ==============================
        lines += """*
* ============ 
* Connectivite 
* ============ 
Zmatrice
"""
        atoms_map = self.connectivity.atoms_map
        lines += "%d\n" % len(atoms_map)
        for i in atoms_map:
            lines += "{} {}\n".format(i, " ".join([str(j) for j in atoms_map[i]]))

        # Atomic charges
        # ==============================
        lines += """*
* =================
* Charges atomiques 
* =================
ModificationChargeDesAtomes e-
"""     
        ChargesMOD = self.ChargesMOD
        lines += "%d\n" % len(ChargesMOD)  # (len(ChargesMOD) + len(self.connectivity.modified_atoms))
        for le in ChargesMOD:
            lines += "%s\n" % le

        # for at in self.connectivity.modified_atoms:
        #     lines += "%d%10.3f\n" % (at, self.connectivity.nodes[at]["charge"])

        # Contribution to dispersion intramolecular
        # ==============================
        lines += """*
* =========================================================
* Contribution de dispersion repulsion en intramoleculaire 
* =========================================================
ContribDispRepIntra
"""     
        ContribDispRepIntra = self.ContribDispRepIntra
        lines += "%d\n" % len(ContribDispRepIntra)
        for le in ContribDispRepIntra:
            lines += "%s\n" % le

        # Save lines
        with open(out, "w") as FATM:
            FATM.write(lines)

        print(f"file name: {out}", end=" - ")
