#!/bin/bash
rm -rf seismograms/rectangle.true.p.mtx
export OMP_NUM_THREADS=1
mpirun -np 4 ./../bin/IFOS "configuration/configuration_true.txt"
