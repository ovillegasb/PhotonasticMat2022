#!/bin/bash

for file in `ls *.sh`
do
    sbatch -A qev@cpu $file
done