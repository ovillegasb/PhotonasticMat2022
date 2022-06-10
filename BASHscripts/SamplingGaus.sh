#!/bin/bash

for i in `ls *.com`
do
    name=${i%.*}
    out=sg16.$name.sh

    rm -vf $out
    touch $out

    echo "#!/bin/bash" >> $out
    echo "#SBATCH --nodes=1 " >> $out
    echo "#SBATCH --ntasks-per-node=1" >> $out
    echo "#SBATCH --cpus-per-task=40" >> $out
    echo "#SBATCH --hint=nomultithread" >> $out
    echo "#SBATCH --job-name=${name}" >> $out
    echo "#SBATCH --input=%x.com" >> $out
    echo "#SBATCH --output=%x.out" >> $out
    echo "#SBATCH --error=%x.o%j" >> $out
    echo "#SBATCH --time=12:00:00" >> $out
    echo "" >> $out
    echo "export jobname=\"${name}\"" >> $out
    echo "" >> $out
    echo "# Manage modules" >> $out
    echo "module purge" >> $out
    echo "module load gaussian/g16-revC01" >> $out
    echo "" >> $out
    echo "export GAUSS_SCRDIR=\$SCRATCH" >> $out
    echo "HW=\`pwd\`" >> $out
    echo "cp \$jobname.com \$SCRATCH " >> $out
    echo "cd \$SCRATCH " >> $out
    echo "g16 -c=\"0-39\" -m=64GB < \$jobname.com" >> $out
    echo "gzip \$jobname.chk" >> $out
    echo "mv \$jobname.chk.gz \$HW" >> $out

done