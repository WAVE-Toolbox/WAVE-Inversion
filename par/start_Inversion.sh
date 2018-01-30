#!/bin/bash
rm -rf seismograms/seismogram.p.mtx
export OMP_NUM_THREADS=1
mpirun -np 4 ./../build/bin/IFOS "configuration/configuration.txt"
