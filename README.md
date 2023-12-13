# Photonastic Materials Repository

Repository dedicated to compile modules and programs created during my PostDoc on photonastic materials.

**Orlando VILLEGAS BELLO**

*2022*

## Install Stamptools

```bash
conda create --name SimMOL scipy numpy pandas
conda activate SimMOL
conda config --add channels conda-forge
conda install --channel conda-forge pymatgen


conda install mdtraj
conda install cython
python setup.py install --prefix="$HOME/.local/"
```


## Modules in repository


| Name      | Code       |
|-----------|------------|
| XYZ2STAMP | python 3.8 |



### Module XYZ2STAMP

Module created to generate input files to STAMP.

    python -m xyz2stamp -option [par]

Force field GAFF (http://ambermd.org/antechamber/gaff.html).