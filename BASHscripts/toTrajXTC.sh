#!/bin/bash

for file in *.xyz
do
	echo $file
	../xyz2gro $file 343 MET 5
done

cat *.gro > traj.gro

gmx trjconv -f traj.gro -o traj.xtc

mv PasDeCalcul_Iteration_00100000.gro confout.gro