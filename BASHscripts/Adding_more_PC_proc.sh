
# Etape I

# N PC + PB
gmx insert-molecules -ci ../../0_files/azoOT.pdb -o azob_box.gro -box 8.0 8.0 8.0 -nmol 2
gmx insert-molecules -ci ../../0_files/btdn.pdb -f azob_box.gro -o box.gro -nmol 2160 -scale 0.4
python ~/.gitproyects/PhotonasticMat/PythonScripts/boxToFAtomes.py -box box.gro -s ../../0_files/{azoOT.pdb,btdn.pdb} -dummyPC


# Aggregate
gmx insert-molecules -ci ../../0_files/azoOT.pdb -o azob_box.gro -box 1.9 1.9 1.9 -nmol 20 -try 1000
gmx editconf -f azob_box.gro -o azob_box_center.gro -bt cubic -c -box 8.0 8.0 8.0
# gmx insert-molecules -ci ../../0_files/btdn.pdb -f azob_box_center.gro -o box.gro -nmol 2160 -scale 0.4
gmx insert-molecules -ci ../../0_files/btdn.pdb -f agg/nvt/confout.gro -o box.gro -nmol 2160 -scale 0.4
python ~/.gitproyects/PhotonasticMat/PythonScripts/boxToFAtomes.py -box box.gro -s ../../0_files/{azoOT.pdb,btdn.pdb} -dummyPC


# Only Polymer
gmx insert-molecules -ci ../0_files/btdn.pdb -o box.gro -nmol 2160 -scale 0.4 -box 8.0 8.0 8.0
python ~/.gitproyects/PhotonasticMat/PythonScripts/boxToFAtomes.py -box box.gro -s ../0_files/btdn.pdb
## add OH



# etape add PC parms and eq NPT
cd iso/2_min
python ~/.gitproyects/PhotonasticMat/PythonScripts/boxToFAtomes.py -f ../1_poly/FATOMES/FAtomes_000500000.in -s ../../../0_files/azoOT.mol2 --addPC --add_OH



### NEW
cd iso/2_min
python ~/.gitproyects/PhotonasticMat/PythonScripts/boxToFAtomes.py -f ../1_poly/FATOMES/FAtomes_000500000.in -s ../../../0_files/azoOT.mol2 --add_OH
cd iso/4_npt
python ~/.gitproyects/PhotonasticMat/PythonScripts/boxToFAtomes.py -f ../1_poly/FATOMES/FAtomes_000500000.in -s ../../../0_files/azoOT.mol2 --addPC