#!/bin/bash
SIMULATIONEXE=~/projects/FDSimulation_LAMA/build/bin/Simulation
MODELEXE=../../build/bin/tools/ModelEttlingerLinie
INVERSIONEXE=../../build/bin/Inversion

# remove old output
rm -f ../seismograms/ettlinger_linie*
rm -f ../gradients/ettlinger_linie*
rm -f ../model/model*
rm -f ../logs/ettlinger_linie*

export SCAI_LOG=ERROR
export SCAI_UNSUPPORTED=IGNORE
export OMP_NUM_THREADS=1

# Create Model
#${MODELEXE} "Input/configuration_true.txt" || { echo 'Model Creation failed' ; exit 1; }

# run WAVE-Simulation to calculate synthetic "field" seismograms
#mpirun -np 12 ${SIMULATIONEXE} "Input/configuration_true.txt" || { echo 'Forward Run failed' ; exit 1; }

# Create Model
${MODELEXE} "Input/configuration.txt" || { echo 'Model Creation failed' ; exit 1; }

# run WAVE-Inversion to invert Model
mpirun -np 4 ${INVERSIONEXE} "Input/configuration.txt"
