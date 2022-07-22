#!/bin/bash

#
# Steps for get OM in local terminal
#

# load modules before
# module load pgi/19.10 gaussian/g16-revC01


inputchk=${1%.*}

# 1
formchk  $jobname.chk

for N in `seq 45 49`
do
    cubegen 1 MO=$N $jobname.fchk $N.cube
done