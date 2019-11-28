#!/bin/bash

# Setup test environment
export OMP_NUM_THREADS=1
export NUM_MPI_PROCESSES=4
export SCAI_UNSUPPORTED=IGNORE
export SCAI_LOG=ERROR

BINDIR="./../build/bin"

INVERSION_EXE="${BINDIR}/Inversion"
UNITTEST_EXE="${BINDIR}/Test_unit"
INTEGRATIONTEST_EXE="${BINDIR}/Test_integration"


# Run unit tests
${UNITTEST_EXE}
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit 1
fi

# Run simulation tests
rm -rf ci/*.ci.*

mpirun -np ${NUM_MPI_PROCESSES} ${INVERSION_EXE} ./ci/configuration_ci.2D.acoustic.txt
 ${INTEGRATIONTEST_EXE} ci/configuration_ci.2D.acoustic.txt
if [ "$?" != "0" ]; then
	echo "Test failed ! "
	exit 2
fi

