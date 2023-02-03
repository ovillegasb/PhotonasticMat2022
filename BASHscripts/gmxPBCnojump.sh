

# Trajetory
# Center asphalten in the box and change the periodicity representation
#echo -e "0 \n" | gmx trjconv -f traj.trr -s run.tpr -o traj_comp.xtc -n index.ndx

echo -e "2 \n 0 \n" | gmx trjconv -f traj_comp.xtc -s run.tpr -o traj_nojump.xtc -center -pbc nojump -n index.ndx

echo -e "2 \n 0 \n" | gmx trjconv -f traj_nojump.xtc -s run.tpr -o traj_nopbc.xtc -center -pbc mol -ur compact -n index.ndx


# Selection from a molecule
# selec a group. The same group for all
echo -e "2 \n 2 \n" | gmx editconf -f confout.gro -o mol.gro -n index.ndx -c
echo -e "2 \n 2 \n" | gmx trjconv -f traj_comp.xtc -s run.tpr -o traj_nojump_mol.xtc -center -pbc nojump -n index.ndx
echo -e "2 \n 2 \n" | gmx trjconv -f traj_nojump_mol.xtc -s run.tpr -o traj_nopbc_mol.xtc -center -pbc mol -ur compact