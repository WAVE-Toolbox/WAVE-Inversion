#!/bin/bash

# Setup test environment
export OMP_NUM_THREADS=1
export NUM_MPI_PROCESSES=4
export SCAI_UNSUPPORTED=IGNORE

SOFI_EXE="./../build/bin/SOFI"
UNITTEST_EXE="./../build/bin/Test_unit"
#INTEGRATIONTEST_EXE="./../build/bin/Test_CompareMisfit"


# Run unit tests
${UNITTEST_EXE}
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit 1
fi

# Run simulation tests
rm -rf ci/*.ci.*

mpirun -np ${NUM_MPI_PROCESSES} ${SOFI_EXE} ci/configuration_ci.2D.acoustic.txt
# ./../bin/Tests/Test_CompareMisfit ci/configuration_ci.2D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	return
fi

