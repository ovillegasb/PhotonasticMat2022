# System preparation
#--------------------

# Locals bin directory: /home/ovillegas/.local/bin
# symbolic links: ln -s path name

# python enviroment
conda activate SimMOL

# STAMP
#---------

# gen XTC trajectory
bash ~/.gitproyects/PhotonasticMat/BASHscripts/gro_to_xtc.sh

# gen trajectory for a molecule
python -m stamptools -d DONNEES.in --mol_traj 0 --traj_type XYZ[or GRO, XTC]

# extract mol properties (Rg, k2, dmax,..)
python -m stamptools -d DONNEES.in --molprop --traj_type XYZ[or GRO, XTC]

# Closest
python -m stamptools -d DONNEES.in -mref 0 --closestDist --reset --traj_type XYZ[or GRO, XTC]

# to convert XYZs files to XTC trajectory gromacs
xyz2gro -d DONNEES.in

# save information about geometry
vmd mol_0_traj.xyz -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_geometry_data.tcl

# RDF from vmd

## XYZ from stamp directly
vmd XYZ/PasDeCalcul_Iteration_0* -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_RDF_data.tcl

## XTC from XYZs STAMP
python ~/GITPROYECTS/PhotonasticMat/PythonScripts/getConfout.py ../DONNEES.in
cat PasDeCalcul__Iteration_*.gro > traj.gro
echo 0 | gmx trjconv -f traj.gro -s confout.gro -o traj_comp.xtc -pbc nojump
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
# Generation confout.gro
python ~/GITPROYECTS/PhotonasticMat/PythonScripts/getConfout.py ../DONNEES.in GRO [None, 0]
# To create a index
echo -e "r 1\n q\n" | gmx make_ndx -f confout.gro -o index.ndx
# Generation mol.gro
echo -e "3\n 3\n" | gmx editconf -f confout.gro -o mol.gro -n index.ndx -c
# new trajectory for a molecule from new gro trajectory
cat PasDeCalcul__Iteration_* >> traj.gro
gmx trjconv -f traj.gro -o traj_comp.xtc

echo -e "3\n 3\n" | gmx trjconv -f traj_comp.xtc -s confout.gro -o traj_comp_mol.xtc -center -pbc nojump -n index.ndx

# Si on dois traslate
# Translate system
echo 0 | gmx trjconv -f traj_comp.xtc -s confout.gro -o traj_translated.xtc -trans -3. -3. 3. -pbc atom



#echo -e "2 \n 2 \n" | gmx trjconv -f traj_comp.xtc -s run.tpr -o traj_nojump_mol.xtc -center -pbc nojump -n index.ndx
#echo -e "2 \n 2 \n" | gmx trjconv -f traj_nojump_mol.xtc -s run.tpr -o traj_nopbc_mol.xtc -center -pbc mol -ur compact




# new trajectory for a molecule from new gro trajectory
cat PasDeCalcul__Iteration_* >> traj.gro
# cp PasDeCalcul__Iteration_10000000.gro confout.gro
python ~/GITPROYECTS/PhotonasticMat/PythonScripts/getConfout.py ../DONNEES.in

gmx trjconv -f traj.gro -o traj_comp.xtc -pbc nojump
echo -e "r 1\nq\n" | gmx make_ndx -f confout.gro -o index.ndx
echo -e "3 \n 3 \n" | gmx editconf -f confout.gro -o mol.gro -n index.ndx -c
echo -e "3\n 3\n" | gmx trjconv -f traj_comp.xtc -s confout.gro -o traj_nojump_mol.xtc -center -pbc nojump -n index.ndx
echo -e "2 \n 2 \n" | gmx trjconv -f traj_nojump_mol.xtc -s mol.gro -o traj_nopbc_mol.xtc -center -pbc mol -ur compact


# UV-VIS analysis
geomSampling_2gaus -f traj_mol_translated.xtc -s mol.gro -b 8500 -e 10500 -dt 1.0 --init_t 2500. -sol PB -isomer trans -N 100
geomSampling_2gaus -s GRO/mol.gro -f GRO/traj_nopbc_mol.xtc -b 500 -e 2500 -sol PB -isomer cis --init_t 500.



# Geometry and RDF DATA in gromacs
vmd confout.gro traj_nopbc.xtc -dispdev text -e ~/GITPROYECTS/PhotonasticMat/VMDscripts/save_geometry_data_GMX.tcl


# MSD
# x label "tau (ps)"
# y label "MSD (nm\\S2\\N)"
gmx msd -s confout.gro -f traj_comp.xtc -n index.ndx -o msd.xvg -xvg none -mol

# Dihedrals
gmx angle -f traj_comp.xtc -n index.ndx -type dihedral -od angdist.xvg -ov angaver.xvg -ot dihtrans.xvg -oh trhisto.xvg -oc dihcorr.xvg


# Translate system
echo 0 | gmx trjconv -f traj_comp.xtc -s confout.gro -o traj_translated.xtc -trans -3. -3. 3. -pbc atom

# New MSD 
# method STAMP
cp `ls PasDeCalcul__Iteration_* | tail -n 1` confout.gro
rm -fv traj.gro
cat PasDeCalcul__Iteration_*.gro > traj.gro
echo 0 | gmx trjconv -f traj.gro -s confout.gro -o traj_comp.xtc -pbc nojump
python ~/.gitproyects/PhotonasticMat/PythonScripts/compute_MSD.py -f traj_comp.xtc -s confout.gro -dt ***0.001*** -o ../msd.csv --select resid 1
rm -v traj.gro

# method GROMACS
rm -fv traj_nojump.gro
echo 0 | gmx trjconv -f traj_comp.xtc -s confout.gro -o traj_nojump.xtc -pbc nojump
python ~/.gitproyects/PhotonasticMat/PythonScripts/compute_MSD.py -f traj_comp.xtc -s confout.gro -dt 0.001 -o ../msd.csv

# MSD from center of mass
# stamp
python ~/.gitproyects/PhotonasticMat/PythonScripts/CMtraj.py STAMP DONNEES.in 0 
echo 0 | gmx trjconv -f traj_cm.gro -s confout_cm.gro -o traj_nojump_cm.xtc -pbc nojump
python ~/.gitproyects/PhotonasticMat/PythonScripts/compute_MSD.py -f traj_nojump_cm.xtc -s confout_cm.gro -dt 1.0 -o msd_cm.csv --select resid 1
# gro
python ~/.gitproyects/PhotonasticMat/PythonScripts/CMtraj.py GRO confout.gro traj_comp.xtc 0
echo 0 | gmx trjconv -f traj_cm.gro -s confout_cm.gro -o traj_nojump_cm.xtc -pbc nojump
python ~/.gitproyects/PhotonasticMat/PythonScripts/compute_MSD.py -f traj_nojump_cm.xtc -s confout_cm.gro -dt 1.0 -o msd_cm.csv --select resid 1


# ADD Oh to FAtomes
mkdir 7_eq3_0
python -m stamptools -f FAtomes.in --add_OH
python -m stamptools -f FAtomes.in --add_H
python -m stamptools -f FAtomes.in --noPBC
# for i in {2..4}; do mkdir 7_eq3_${i}; cd 7_eq3_${i}/; cp ../6_prod_${i}/FAtomes.in .; python -m stamptools -f FAtomes.in --add_OH; cp ../7_eq3_0/DONNEES.in .; sbatch ~/scripts/stamp.sh DONNEES.in; cd ..; done

mkdir 8_solOH_NPT_0
cd 8_solOH_NPT_0
cp ../7_eq3_0/FAtomes_modif.in FAtomes.in  # Modificar reatividad
cp ../8_solOH_NPT_0/DONNEES.in .
cp ../8_solOH_NPT_0/stamp.sh .
ls ../7_eq3_0/PROTS/* | tail -n 1 > DerniereProtection

cat Stamp.log | grep "Minimum global"


# change timestep
gmx trjconv -f traj_nopbc_mol.xtc -s mol.gro -o traj_nopbc_mol_1fs.xtc -dt 1.0

# Gen sampling from a system
python -m stamptools -d DONNEES.in --mol_d_aa -b 1500 -e 2000 -i 10 -mref 0 --rcutoff 0.12 --out_folder sampling

# UV-VIS analysis
geomSampling_2gaus -s mol.gro -f traj_mol_translated.xtc -b 8500 -e 10500 -dt 1.0 --init_t 2500. -sol PB -isomer trans -N 100
geomSampling_2gaus -s GRO/mol.gro -f GRO/traj_nopbc_mol.xtc -b 500 -e 2500 -sol PB -isomer cis --init_t 500.

# For N PC in the system
bash ~/.gitproyects/PhotonasticMat/BASHscripts/gro_to_xtc.sh
for i in {0..N-1}; do python -m stamptools -d DONNEES.in --mol_traj $i --traj_type XTC; done
for i in {0..N-1}; do vmd mol_traj_$i.xyz -dispdev text -e ~/.gitproyects/PhotonasticMat/VMDscripts/save_geometry_N_pc.tcl; done

for i in {0..1}; do python -m stamptools -d DONNEES.in --mol_traj $i --traj_type GRO; vmd mol_traj_$i.xyz -dispdev text -e ~/.gitproyects/PhotonasticMat/VMDscripts/save_geometry_N_pc.tcl; done

vmd XTC/confout.gro XTC/traj_comp.xtc -dispdev text -e ~/.gitproyects/PhotonasticMat/VMDscripts/save_hbonds.tcl

# Aggregates
python -m stamptools -d DONNEES.in --aggregate 0:19 --traj_type GRO