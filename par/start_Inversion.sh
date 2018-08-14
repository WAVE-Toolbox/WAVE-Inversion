#!/bin/bash
rm -rf seismograms/seismogram.p.mtx
rm -f gradients/*.mtx
rm -f model/*It*
export OMP_NUM_THREADS=1
export SCAI_UNSUPPORTED=IGNORE
mpirun -np 4 ./../build/bin/IFOS "configuration/configuration.txt"
