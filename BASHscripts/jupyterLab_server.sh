#!/bin/bash

source /home/${USER}/.bashrc
conda -V
conda activate SimMOL
echo "CONDA ENV: ${CONDA_DEFAULT_ENV}"
echo "date: `date`"

# Init jupyter-lab
#TOKEN=$(jupyter lab --no-browser --port=8888 | awk '/http/{print $NF}')

#echo "Token generado: $TOKEN"

# Send process to backgroud
# nohup jupyter lab --no-browser --port=8888 --ip=0.0.0.0 &> /dev/null &
# Normal: 8888 ---> 8898
nohup jupyter lab --no-browser --port=8898 --ip=0.0.0.0 &> /tmp/nohup.out & 
#nohup jupyter lab --no-browser --port=8888 &> /tmp/nohup.out &


#cat nohup.out

#tail -n 2 nohup.out | head -n 1 | sed 's/yoda:8888/localhost:8889/g'

echo "Token using:"
echo "tail -n 2 /tmp/nohup.out | head -n 1 | sed 's/yoda:8898/localhost:8889/g'"

