#!/bin/bash
SIMULATIONEXE=${GPI_ROOT}/bin/Simulation
MODELEXE=../../build/bin/tools/ToyModel
IFOSEXE=../../build/bin/IFOS

# remove old output
rm -f ../seismograms/ToyExample*
rm -f ../gradients/ToyExample*
rm -f ../model/ToyModel3D*
rm -f ../logs/ToyExample*

#export SCAI_LOG=ERROR
export SCAI_UNSUPPORTED=IGNORE
export OMP_NUM_THREADS=1

# Create Model
${MODELEXE} "Input/configuration_true.txt" || { echo 'Model Creation failed' ; exit 1; }

# run WAVE-Simulation to calculate synthetic "field" seismograms
mpirun -np 8 ${SIMULATIONEXE} "Input/configuration_true.txt" || { echo 'Forward Run failed' ; exit 1; }

# run IFOS to invert Model
mpirun -np 8 ${IFOSEXE} "Input/configuration.txt"

echo 'finished'
