#!/bin/bash

source /home/${USER}/.bashrc
conda -V
conda activate SimMOL
echo "CONDA ENV: ${CONDA_DEFAULT_ENV}"
echo "date: `date`"

# Init jupyter-lab
TOKEN=$(jupyter lab --no-browser --port=8888 | awk '/http/{print $NF}')

echo "Token generado: $TOKEN"

# Send process to backgroud
nohup jupyter lab --no-browser &> /dev/null &

