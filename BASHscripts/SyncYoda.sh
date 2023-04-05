#!/bin/bash

rsync -avrhP \
--exclude 'exp.azob.stamp' \
--exclude 'test' \
--exclude=".*" \
--exclude="#*" \
--include=".bashrc" \
--delete \
ovillegas@10.32.153.46:/home/ovillegas /run/media/ovillegas/My\ Passport/YODA/


# --dry-run \