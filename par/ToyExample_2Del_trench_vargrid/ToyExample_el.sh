#!/bin/bash
SOFIEXE=${GPI_ROOT}/bin/SOFI
MODELEXE=../../build/bin/tools/ModelEttlingerLinie
IFOSEXE=../../build/bin/IFOS

# remove old output
rm -f ../seismograms/ettlinger_linie*
rm -f ../gradients/ettlinger_linie*
rm -f ../model/model*
rm -f ../logs/ettlinger_linie*

export SCAI_LOG=ERROR
export SCAI_UNSUPPORTED=IGNORE
export OMP_NUM_THREADS=1

# Create Model
${MODELEXE} "Input/configuration_true.txt" || { echo 'Model Creation failed' ; exit 1; }

# run SOFI to calculate synthetic "field" seismograms
mpirun -np 4 ${SOFIEXE} "Input/configuration_true.txt" || { echo 'Forward Run failed' ; exit 1; }

# Create Model
${MODELEXE} "Input/configuration.txt" || { echo 'Model Creation failed' ; exit 1; }

# run IFOS to invert Model
mpirun -np 4 ${IFOSEXE} "Input/configuration.txt"
