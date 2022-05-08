#!/bin/bash

conda create --name SimMOL scipy jupyterlab numpy matplotlib seaborn pandas
conda activate SimMOL
conda install --channel conda-forge pymatgen
