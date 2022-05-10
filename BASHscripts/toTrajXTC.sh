#!/bin/bash

conf=$1
NMOL=$2
RES=$3
NAT=$4

for file in *.xyz
do
	echo $file
	xyz2gro $file $NMOL $RES $NAT
done

cat *.gro > traj.gro

gmx trjconv -f traj.gro -o traj.xtc

mv "${conf%.*}.gro" confout.gro

# rm *.xyz
