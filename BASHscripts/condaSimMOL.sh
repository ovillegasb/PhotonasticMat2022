#!/bin/bash

conda create --name SimMOL scipy jupyterlab numpy matplotlib seaborn pandas
conda activate SimMOL
conda install --channel conda-forge pymatgen
conda install nglview -c conda-forge

conda config --add channels conda-forge
conda install mdanalysis
conda install -c anaconda networkx

# Others
sudo dnf install pandoc
conda install ipympl  # %matplotlib widget
conda install -c conda-forge nodejs
conda install -c conda-forge/label/gcc7 nodejs
conda install -c conda-forge/label/cf201901 nodejs
conda install -c conda-forge/label/cf202003 nodejs
conda install -n base -c conda-forge jupyterlab_widgets
conda install -n SimMOL -c conda-forge ipywidgets
conda install -c conda-forge mpi4py
conda install memory_profiler
conda install statsmodels
pip install pyppeteer
conda install -c conda-forge mdtraj
# conda install -c auto threadpool
# conda install -c auto multiprocessing


# Enviroment AmberTools
conda create --name AmberTools22
conda activate AmberTools22
conda install -c conda-forge ambertools=22 compilers
conda install -c conda-forge acpype
# acpype -i test/mol2/hpt.mol2 -a gaff2 -o gmx


# To convert Jupyter to slides
jupyter nbconvert --to slides --no-input --post serve presentation.slides.ipynb 