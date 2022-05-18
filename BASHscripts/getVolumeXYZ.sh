#!/bin/bash


rm -vf box.dat
touch box.dat


for file in *.xyz
do
    # get second line from xyz file
    head -n 2 $file | sed -n 2p >> box.dat
done