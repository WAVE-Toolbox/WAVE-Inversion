#!/bin/bash

# Setup test environment
export OMP_NUM_THREADS=1
export NUM_MPI_PROCESSES=4
export SCAI_UNSUPPORTED=IGNORE

# Run unit tests
./../bin/Tests/Test_unit
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit 1
fi

# Run simulation tests
rm -rf ci/*.ci.*

mpirun -np ${NUM_MPI_PROCESSES} ./../bin/IFOS ci/configuration_ci.2D.acoustic.txt
 ./../bin/Tests/Test_CompareMisfit ci/configuration_ci.2D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

