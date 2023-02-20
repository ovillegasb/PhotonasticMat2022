# System preparation
#--------------------

# Locals bin directory: /home/ovillegas/.local/bin
# symbolic links: ln -s path name

# python enviroment
conda activate SimMOL

# STAMP
#---------

# gen trajectory for a molecule
python -m stamptools -d DONNEES.in --mol_traj 0

# to convert XYZs files to XTC trajectory gromacs
xyz2gro -d DONNEES.in

# save information about geometry
vmd mol_0_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_geometry_data.tcl

# RDF from vmd

## XYZ from stamp directly
vmd XYZ/PasDeCalcul_Iteration_0* -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_RDF_data.tcl

## XTC from XYZs STAMP
vmd GRO_files/confout.gro GRO_files/traj_comp.xtc -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_RDF_data_GMX2.tcl

# Geometry analysis

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

# Geometry and RDF DATA in gromacs
vmd confout.gro traj_nopbc.xtc -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_geometry_data_GMX.tcl