CC=g++
CFLAGS=-std=c++11 -fopenmp -Wall -Wextra -Werror 
IFLAGS=-I${GPI_ROOT}/src/ -I${SCAI_ROOT}/include  -DSCAI_ASSERT_LEVEL_DEBUG -DSCAI_LOG_LEVEL_DEBUG -DSCAI_TRACE_OFF 
#-I../SOFI/src
LDFLAGS=-L${GPI_ROOT}/bin -lSOFI -L${SCAI_ROOT}/lib  -lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common 
#LDFLAGS= -L${SCAI_ROOT}/lib -lscai_lama -lscai_solver -lscai_dmemo -lscai_sparsekernel -lscai_utilskernel -lscai_blaskernel -lscai_kregistry -lscai_hmemo -lscai_tasking -lscai_tracing -lscai_logging -lscai_common -L ~/projects/FDSimulation_LAMA/bin -lSOFI

# Flags for googletest framework used for unit tests
IFLAGSTEST = -I${GTEST_DIR}/include
LDFLAGSTEST = -L${GTEST_DIR}/ -lgtest -lgtest_main
CFLAGSTEST = -isystem -pthread 

######################
# Install directory
######################
BINDIR=../bin

LIB_TARGET=$(BINDIR)/libIFOS.a
IFOS_TARGET=$(BINDIR)/IFOS
IFOS_UNITTEST_TARGET=$(BINDIR)/Tests/Test_unit
#IFOS_INTEGRATIONTEST_TARGET=$(BINDIR)/Tests/Test_CompareSeismogram

################################################################
################################################################
# No input parameters are required below this line
################################################################
################################################################

######################
# General targets
######################
all: ifos tests

install: clean ifos 

tests: unitTest
#tests: integrationTest unitTest

lib: $(LIB_TARGET)

ifos: $(IFOS_TARGET)

unitTest: $(IFOS_UNITTEST_TARGET)

#integrationTest: $(IFOS_INTEGRATIONTEST_TARGET)

.PHONY : clean
clean: clean-bin clean-obj

######################
# Pattern rule
######################
.cpp.o:
	$(CC) $(CFLAGS) -c  $< -o $@ $(IFLAGS)

SUBDIRS := $(wildcard */)

######################
# Locate CPP files for compilation
######################
IFOS_LIB_CPP := $(wildcard $(addsuffix */*.cpp,$(SUBDIRS))) $(wildcard **/*.cpp)
IFOS_UNITTEST_CPP := $(wildcard $(addsuffix */*UnitTest.cpp,$(SUBDIRS))) $(wildcard */*UnitTest.cpp)
#IFOS_INTEGRATIONTEST_CPP := $(wildcard **/*Test_CompareSeismogram.cpp)

######################
# Static IFOS library 
######################
IFOS_LIB_SCR := $(filter-out $(IFOS_UNITTEST_CPP) $(IFOS_INTEGRATIONTEST_CPP), $(IFOS_LIB_CPP))
IFOS_LIB_OBJ := $(IFOS_LIB_SCR:%.cpp=%.o)

$(LIB_TARGET): $(IFOS_LIB_OBJ)
	ar rcuvs $@ $^
	ranlib $@

#####################
# IFOS main program
#####################
IFOS_MAIN_CPP := IFOS.cpp
IFOS_MAIN_OBJ := $(IFOS_MAIN_CPP:%.cpp=%.o)

$(IFOS_TARGET): $(IFOS_MAIN_OBJ) $(LIB_TARGET)
	$(CC) $(CFLAGS) $^ -o $@ $(IFLAGS) $(LDFLAGS)

######################
# Unit tests
######################
IFOS_UNITTEST_OBJ = $(IFOS_UNITTEST_CPP:%.cpp=%.o)
$(IFOS_UNITTEST_OBJ): CFLAGS+=$(CFLAGSTEST)
$(IFOS_UNITTEST_OBJ): LDFLAGS +=$(LDFLAGSTEST)
$(IFOS_UNITTEST_OBJ): IFLAGS +=$(IFLAGSTEST)

$(IFOS_UNITTEST_TARGET): $(IFOS_UNITTEST_OBJ) $(LIB_TARGET)
	mkdir -p $(BINDIR)/Tests/
	$(CC) $(CFLAGS) $(CFLAGSTEST) $^ -o $@ $(IFLAGS) $(IFLAGSTEST) $(LDFLAGS) $(LDFLAGSTEST)

######################
# Integration tests
######################
#IFOS_INTEGRATIONTEST_OBJ := $(IFOS_INTEGRATIONTEST_CPP:%.cpp=%.o)

#$(IFOS_INTEGRATIONTEST_TARGET): $(IFOS_INTEGRATIONTEST_OBJ) $(LIB_TARGET)
#	mkdir -p $(BINDIR)/Tests/
#	$(CC) $(CFLAGS) $^ -o $@ $(IFLAGS) $(LDFLAGS)

######################
# Cleaning
######################
.PHONY : clean-bin
clean-bin:
	rm -rf $(BINDIR)/IFOS
	rm -rf $(BINDIR)/libIFOS.a


.PHONY : clean-obj
clean-obj:
	rm -rf $(IFOS_LIB_OBJ) $(IFOS_MAIN_OBJ)

