# Photonastic Materials Repository

Repository dedicated to compile modules and programs created during my PostDoc on photonastic materials.

**Orlando VILLEGAS BELLO**

*2022*

## Install Stamptools

```bash
conda create --name SimMOL scipy jupyterlab numpy matplotlib seaborn pandas
conda activate SimMOL
conda config --add channels conda-forge
conda install --channel conda-forge pymatgen

conda install pymatgen
conda install mdtraj
conda install mdanalysis
conda install networkx
conda install statsmodels
python setup.py install

# Install MolCraft
cd molcraft/
python setup.py install

# Install StampTools
cd PhotonasticMat/
python setup.py install
```

## Scripts included in `PhotonasticMat` repertoire

### `./BASHscripts/`

+   `Adding_more_PC_proc.sh`. Procedure to add more than one PC to a STAMP system.
+   `AmberToolsEnv.sh`. Configuring the AmberTools environment.
+   `build_GROMACS.sh`. Cmake instruction to compile gromacs.
+   `condaSimMOL.sh`. Useful commands used to configure with SimMOL.
+   `SyncYoda.sh` and `SyncBettyBoop.sh`. Files to synchronize the server with a drive such as a removable hard disk.
+   `jupyterLab_server.sh`. Run jupyterlab in the server backgound.
+   `gro_to_xtc.sh`. Procedure to transform gro files to xtc.


### `./GMXfiles/`

It contains the gaff force field in gromacs format and mdp calculation configuration files.

### `./VMDscripts/`

Contains routines and functions to be used in VMD.


### `./PythonScripts/`

It contains files and routines written in Python.

+   `xyz2gro.py `. Script used to convert STAMP xyz trajectory to GROMACS gro trajectory.
+   `getConfout.py`. Generate a confout.gro gromacs file removing the periodical conditions.
+   `geomSampling_2gaus.py`. Generates a series of gaussian input files sampling a trajectory.
+   `genGaussFiles.py`. Generate Gaussian inputs using ONIOM methods.
+   `CMtraj.py`. It generates a trajectory of centers of masses.
+   `compute_MSD.py`. Calculates the MSD of a system.