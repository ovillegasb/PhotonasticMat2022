#!/bin/bash

# Script used to calculate the partial charges of a
# molecule from a .pdb file using antechamber.

# conda activate AmberTools22

antechamber -i ${1} -fi pdb -o ${1%.*}.mol2 -fo mol2 -c bcc -nc 0