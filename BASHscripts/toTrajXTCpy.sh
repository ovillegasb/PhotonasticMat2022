#!/bin/bash

# final frame
conf=$1
# Number of residues
NMOL=$2
# Residue name
RES=$3
# Number of atoms per residue
NAT=$4
# Piston velocity (m/s) positive
VEL=$5

# echo "${NMOL}-${RES}-${NAT}"
# echo $VEL

for file in *.xyz
do
	echo $file
	python ~/GITPROYECTS/PhotonasticMat/PythonScripts/xyz2gro.py  $file -r $NMOL-$RES-$NAT -p $VEL
done

cat *.gro > traj.gro

gmx trjconv -f traj.gro -o traj.xtc

mv "${conf%.*}.gro" confout.gro

# rm *.xyz
