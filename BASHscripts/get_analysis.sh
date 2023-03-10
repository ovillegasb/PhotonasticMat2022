#!/bin/bash

# conda activate SimMOL

if [[ ${CONDA_DEFAULT_ENV}!="SimMOL" ]]
then
    source /home/${USER}/.bashrc
    conda activate SimMOL
fi

echo "CONDA ENV: ${CONDA_DEFAULT_ENV}"

for isomer in cis trans
do
    cd $isomer/
    for r in 0 1 2 3 4
    do
        cd 6_prod_$r/
        
        # UV-VIS sampling analysis
        ## python -m stamptools -d DONNEES.in --mol_traj 0 --format gro
        ## geomSampling_2gaus -s mol.gro -f traj_comp_mol_0.xtc -b 500 -e 2500 -sol PB -isomer $isomer --init_t 500.

        # Geometry analysis
        python -m stamptools -d DONNEES.in --mol_traj 0
        vmd mol_0_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_geometry.tcl

        # RDF all atoms analysis
        ##vmd GRO/*.gro -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_RDF.tcl

        # Polymer properties
        # python -m stamptools -d DONNEES.in --molprop
        # python -m stamptools -d DONNEES.in -mref 0 --closestDist

        cd ../
    done
    cd ../
done
