#!/bin/bash

conda create --name SimMOL scipy jupyterlab numpy matplotlib seaborn pandas
conda activate SimMOL
conda install --channel conda-forge pymatgen
conda install nglview -c conda-forge

conda config --add channels conda-forge
conda install mdanalysis

# Others
sudo dnf install pandoc
conda install -c conda-forge nodejs
conda install -c conda-forge/label/gcc7 nodejs
conda install -c conda-forge/label/cf201901 nodejs
conda install -c conda-forge/label/cf202003 nodejs
conda install -n base -c conda-forge jupyterlab_widgets
conda install -n SimMOL -c conda-forge ipywidgets
