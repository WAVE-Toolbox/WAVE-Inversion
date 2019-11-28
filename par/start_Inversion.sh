#!/bin/bash
BINDIR="./../build/bin"

export OMP_NUM_THREADS=1
export SCAI_UNSUPPORTED=IGNORE

mpirun -np 4 "${BINDIR}/Inversion" "configuration/configuration.txt"
