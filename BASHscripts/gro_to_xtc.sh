
source /home/${USER}/.bashrc
conda activate SimMOL

# In folder system

if [ ! -d XTC ]
then
    mkdir XTC
else
    rm -vf XTC/*
fi

cd XTC

# gen confout.gro
python ~/.gitproyects/PhotonasticMat/PythonScripts/getConfout.py ../DONNEES.in GRO None   #   XYZ[or GRO] resid_pc

# trajectory
if [ -f traj.gro ]
then
    rm -v traj.gro
fi

cat ../GRO/PasDeCalcul_* > traj.gro
echo 0 | gmx trjconv -f traj.gro -s confout.gro -o traj_comp.xtc
rm -v traj.gro
 
cd ..
echo "Finished!"
