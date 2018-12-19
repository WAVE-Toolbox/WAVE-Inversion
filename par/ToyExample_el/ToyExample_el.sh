#!/bin/bash
SOFIEXE=${GPI_ROOT}/bin/SOFI
MODELEXE=../../build/bin/tools/ToyModel
IFOSEXE=../../build/bin/IFOS

# remove old output
rm -f ../seismograms/ToyExample*
rm -f ../gradients/ToyExample*
rm -f ../model/ToyModel*
rm -f ../logs/ToyExample*

export SCAI_LOG=ERROR
export SCAI_UNSUPPORTED=IGNORE
export OMP_NUM_THREADS=1

# Create Model
${MODELEXE} "Input/configuration_true.txt" || { echo 'Model Creation failed' ; exit 1; }

# run SOFI to calculate synthetic "field" seismograms
mpirun -np 4 ${SOFIEXE} "Input/configuration_true.txt" || { echo 'Forward Run failed' ; exit 1; }

# run IFOS to invert Model
mpirun -np 4 ${IFOSEXE} "Input/configuration.txt"
