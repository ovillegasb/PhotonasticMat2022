#!/bin/bash

# conda activate SimMOL

if [[ ${CONDA_DEFAULT_ENV}!="SimMOL" ]]
then
    source /home/${USER}/.bashrc
    conda activate SimMOL
fi

echo "CONDA ENV: ${CONDA_DEFAULT_ENV}"
echo "============================="
sleep 5

# COLORS
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'


# Somes functions

make_traj_gmx () {
    echo "Trajectory gromacs"
    sleep 5
    cp $(ls PasDeCalcul__Iteration_* | tail -n 1) confout.gro
    cat PasDeCalcul__Iteration_* > traj.gro
    gmx trjconv -f traj.gro -o traj_comp.xtc
    echo -e "r 1\n name 3 PHO\n q\n" | gmx make_ndx -f confout.gro -o index.ndx
}

start=$(date +%s)

for iso in cis trans
do
    cd $iso/
    for r in 0 1 2 3 4
    do
        cd 6_prod_$r/
        echo "Isomer: ${iso} - replica ${r}"
        
        # UV-VIS sampling analysis
	echo -e "${BLUE}UV-VIS sampling analysis${NC}"
	sleep 5
        python -m stamptools -d DONNEES.in --mol_traj 0 --format gro
        /home/ovillegas/.local/bin/geomSampling_2gaus -s mol.gro -f traj_comp_mol_0.xtc -dt 1.0  -b 500 -e 2500 -sol PB -isomer $iso --init_t 500.
	sleep 5

        # Geometry analysis
	echo -e "${BLUE}Geometry analysis${NC}"
	sleep 5
	python -m stamptools -d DONNEES.in --mol_traj 0
        vmd mol_0_traj.xyz -dispdev text -e ~/.gitproyects/PhotonasticMat/VMDscripts/save_geometry.tcl
	sleep 5

        # RDF all atoms analysis
	echo -e "${BLUE}RDF all atoms analysis${NC}"
	sleep 5
        vmd GRO/PasDeCalcul__Iteration_*.gro -dispdev text -e ~/.gitproyects/PhotonasticMat/VMDscripts/save_RDF.tcl
	sleep 5

        # Polymer properties
        # python -m stamptools -d DONNEES.in --molprop
        # python -m stamptools -d DONNEES.in -mref 0 --closestDist

        # GRO files
        cd GRO
        make_traj_gmx
        ## MSD
	echo -e "${BLUE}MSD analysis${NC}"
	sleep 5
        echo "PHO" | gmx msd -s confout.gro -f traj_comp.xtc -n index.ndx -o ../msd.xvg -xvg none -mol
        cd ..
	sleep 5

	sbatch ~/scripts/spython.sh python -m stamptools -d DONNEES.in --molprop
	sbatch ~/scripts/spython.sh python -m stamptools -d DONNEES.in -mref 0 --closestDist
        echo "============================="
        cd ../
    done
    cd ../
done
end=$(date +%s)

echo "============================="
echo "Elapsed Time: $(($end-$start)) seconds"
echo -e "${RED}Finish!${NC}"
