#!/bin/bash
myid=${OMPI_COMM_WORLD_RANK}

mkdir myLJMD$myid
cd myLJMD$myid

cp ../01_baseline/ljmd-f.x .
cp ../00_input/argon_108.* .

./ljmd-f.x < argon_108.inp | tee output.txt

