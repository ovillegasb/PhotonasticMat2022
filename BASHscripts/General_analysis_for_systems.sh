# System preparation
#--------------------

# Locals bin directory: /home/ovillegas/.local/bin
# symbolic links: ln -s path name

# python enviroment
conda activate SimMOL

# STAMP



# GROMACS
#---------

# preaparation:

# trajectory for all molecules
echo -e "q\n" | gmx make_ndx -f confout.gro -o index.ndx
echo -e "2 \n 0 \n" | gmx trjconv -f traj_comp.xtc -s run.tpr -o traj_nojump.xtc -center -pbc nojump -n index.ndx
echo -e "2 \n 0 \n" | gmx trjconv -f traj_nojump.xtc -s run.tpr -o traj_nopbc.xtc -center -pbc mol -ur compact -n index.ndx

# trajectory for a molecule
echo -e "2 \n 2 \n" | gmx editconf -f confout.gro -o mol.gro -n index.ndx -c
echo -e "2 \n 2 \n" | gmx trjconv -f traj_comp.xtc -s run.tpr -o traj_nojump_mol.xtc -center -pbc nojump -n index.ndx
echo -e "2 \n 2 \n" | gmx trjconv -f traj_nojump_mol.xtc -s run.tpr -o traj_nopbc_mol.xtc -center -pbc mol -ur compact

# UV-VIS analysis
geomSampling_2gaus -f traj_nopbc_mol.xtc -s mol.gro -b 0 -e 2500 -r 0 -sol CHX -isomer cis

# Geometry DATA in gromacs
vmd confout.gro traj_nopbc.xtc -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_geometry_data_GMX.tcl