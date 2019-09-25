###########################################################################
#
#     Start the build of IFOS++ by the following command:
# 
#     SCAI_DIR=<path_to_lama> GTEST_ROOT=<path_to_googletest> cmake ../src 
#
#     make 
#     make install
#
#     default installation directory is the build directory
#    
###########################################################################

project( IFOS )

cmake_minimum_required ( VERSION 3.0 )

if ( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif ()


## Set default installation directory same as build directory

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
    set ( CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )


### get git commit ###
execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
add_definitions("-DBUILD_COMMIT=${GIT_COMMIT}")


####################################################
#  Find SOFI package (mandatory)                   #
####################################################

set ( SOFI_PATH_HINTS "${SOFI_DIR}" "${SOFI_HOME}" "$ENV{SOFI_HOME}" )

if ( ${CMAKE_VERSION} VERSION_GREATER 3.11)
    ## use new policy to configure with SCAI_ROOT, suppress warning
    cmake_policy( SET CMP0074 NEW )
else ()
    set ( SOFI_PATH_HINTS ${SOFI_PATH_HINTS} "${SOFI_ROOT}" "$ENV{SOFI_ROOT}" )
endif ()

find_package( SOFI REQUIRED CONFIG 
              NO_CMAKE_PATH 
              NO_SYSTEM_ENVIRONMENT_PATH 
              PATHS ${SOFI_PATH_HINTS} )

## take over compilfer flags as used by SOFI

set( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SOFI_CXX_FLAGS} )

## provide SOFI include directories + library for IFOS

include_directories( ${SOFI_INCLUDE_DIRS} )
set( IFOS_used_libs ${SOFI_LIBRARIES} )

####################################################
#  Find SCAI/LAMA package (mandatory)              #
####################################################

# take SOFI_SCAI_ROOT as additional hint so same installation is used

set ( SCAI_PATH_HINTS "${SOFI_SCAI_ROOT}" )


find_package( SCAI REQUIRED CONFIG 
              NO_CMAKE_PATH 
              NO_SYSTEM_ENVIRONMENT_PATH 
              PATHS ${SCAI_PATH_HINTS} )

## provide SCAI/LAMA include directories + library for IFOS

include_directories( ${SCAI_INCLUDE_DIRS} )
set( IFOS_used_libs ${IFOS_used_libs} ${SCAI_LIBRARIES} )

####################################################
#  Find Geographer package (optional)              #
####################################################

# take SOFI_Geographer_ROOT as additional hint so same installation is used

set ( Geographer_PATH_HINTS 
        "${Geographer_DIR}" 
        "${Geographer_HOME}"
        "$ENV{Geograhper_HOME}"
        "${SOFI_Geographer_ROOT}" )

if ( ${CMAKE_VERSION} VERSION_GREATER 3.11)
else()
    set ( Geographer_PATH_HINTS ${Geographer_PATH_HINTS} "${Geographer_ROOT}" "$ENV{Geographer_ROOT}" )
endif()

find_package( Geographer CONFIG
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH 
              PATHS ${Geographer_PATH_HINTS} )

if ( ${Geographer_FOUND} )
 
    ## enalble conditional compilation of code using Geographer
    add_definitions( "-DUSE_GEOGRAPHER" )

    ## provide Geographer include directories + library for IFOS
    include_directories( ${Geographer_INCLUDE_DIRS} )
    set( IFOS_used_libs ${Geographer_LIBRARIES} ${IFOS_used_libs} )

endif ()

###################################################
#  Collect source files of IFOS library            #
#  and copy header files to install                #
####################################################

set ( IFOS_LIB_DIRS 
	Common
	Gradient
	Misfit
	Optimization
	Preconditioning
	SourceEstimation
	StepLengthSearch
	Taper
	Workflow
	ZeroLagCrossCorrelation
         )

foreach ( dir ${IFOS_LIB_DIRS} ) 

    file( GLOB_RECURSE files ${dir}/*.cpp ${dir}/*.hpp )
    set ( IFOS_FILES ${IFOS_FILES} ${files} )

    file( GLOB_RECURSE files RELATIVE ${CMAKE_SOURCE_DIR}/${dir} ${dir}/*.hpp )
    foreach ( file ${files} )
	    configure_file( ${dir}/${file} include/${dir}/${file} COPYONLY )
    endforeach ()

endforeach ()

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)


####################################################
#  IFOS library                                    #
####################################################

include_directories( Common Gradient Misfit Optimization Preconditioning SourceEstimation StepLengthSearch Taper Workflow ZeroLagCrossCorrelation )

add_library( IFOS ${IFOS_FILES} )

target_link_libraries( IFOS ${IFOS_used_libs} )

install( TARGETS IFOS DESTINATION lib )

####################################################
#  IFOS executable                                 #
####################################################

add_executable( prog IFOS.cpp )

target_link_libraries( prog IFOS ${IFOS_used_libs} )

set_target_properties( prog PROPERTIES OUTPUT_NAME IFOS )

install( TARGETS prog DESTINATION bin )

####################################################
#  Unit Test                                       #
####################################################

find_package( GTest )

if ( ${GTEST_FOUND} )

  file( GLOB IFOS_UNIT_TEST_SOURCES Tests/UnitTest/*.cpp )

  add_executable( utest ${IFOS_UNIT_TEST_SOURCES} )

  target_include_directories( utest PRIVATE Common Gradient Misfit Optimization Preconditioning SourceEstimation StepLengthSearch Taper Workflow ZeroLagCrossCorrelation  ${GTEST_INCLUDE_DIR} )
  target_link_libraries( utest ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} IFOS ${IFOS_used_libs} )

  set_target_properties( utest PROPERTIES OUTPUT_NAME Test_unit )

  install( TARGETS utest DESTINATION bin )

else()

   message( WARNING "Skip configuration for utest, please set GTEST_ROOT to googletest installation directory" )

   add_custom_target (
       utest
       COMMAND echo "ATTENTION: cannot build utest as googletest not found, set GTEST_ROOT" 
   )

endif ()


#  Integration Test                                #
####################################################

set ( IFOS_INTEGRATION_TEST_SOURCES
   Tests/IntegrationTest/Test_CompareMisfit.cpp
)

add_executable( itest ${IFOS_INTEGRATION_TEST_SOURCES} )

target_link_libraries( itest IFOS ${IFOS_used_libs} )

set_target_properties( itest PROPERTIES OUTPUT_NAME Test_integration )

install( TARGETS itest DESTINATION bin )

#####################################################
##  Create Model                                    #
#####################################################

ADD_CUSTOM_TARGET(model)
ADD_DEPENDENCIES(model TwoLayer)

add_executable( TwoLayer Tools/CreateModel/TwoLayer.cpp )

target_link_libraries( TwoLayer IFOS ${IFOS_used_libs} )

set_target_properties( TwoLayer PROPERTIES OUTPUT_NAME TwoLayer )

install( TARGETS TwoLayer DESTINATION bin/tools )

add_executable( ToyModel Tools/CreateModel/ToyModel.cpp )

target_link_libraries( ToyModel IFOS ${IFOS_used_libs} )

set_target_properties( ToyModel PROPERTIES OUTPUT_NAME ToyModel )

install( TARGETS ToyModel DESTINATION bin/tools )


#####################################################
##  Doxygen documentation                           #
#####################################################

find_package( Doxygen )

message ( STATUS "DOXYGEN_FOUND=${DOXYGEN_FOUND}, DOXYGEN_EXECUTABLE=${DOXYGEN_EXECUTABLE}" )

if ( DOXYGEN_FOUND )

    message( [status] "${CMAKE_SOURCE_DIR}/../doc/doxygen/Doxyfile")

    set ( DOXYGEN_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen" )

    configure_file ( "${CMAKE_SOURCE_DIR}/../doc/doxygen/Doxyfile" "${DOXYGEN_BINARY_DIR}/Doxyfile" @ONLY )

    file ( MAKE_DIRECTORY ${DOXYGEN_BINARY_DIR}/html )

    add_custom_command (
        OUTPUT ${DOXYGEN_BINARY_DIR}/html/index.html
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_BINARY_DIR}/Doxyfile
        DEPENDS ${DOXYGEN_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${DOXYGEN_BINARY_DIR}
    )

    add_custom_target (
        doc
        DEPENDS ${DOXYGEN_BINARY_DIR}/html/index.html
        COMMENT "Creating doxygen doc."
    )

    install ( DIRECTORY   ${DOXYGEN_BINARY_DIR}/html 
              DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc )

else ()

    add_custom_target (
        doc
        COMMAND echo "ATTENTION: doxygen not found, cannot build documentation" 
    )
 
endif ()

####################################################
#  Guide (latex/pdf) document                      #
####################################################

find_package ( LATEX COMPONENTS PDFLATEX BIBTEX )

message ( STATUS "LATEX_FOUND=${LATEX_FOUND}, PDFLATEX_COMPILER=${PDFLATEX_COMPILER}, BIBTEX_COMPILER=${BIBTEX_COMPILER}" )

if ( LATEX_FOUND )

   set ( GUIDE_BINARY_DIR "${CMAKE_BINARY_DIR}/doc/guide" )
   configure_file ( "${CMAKE_SOURCE_DIR}/../doc/guide/WAVE_guide_EN.tex" "${GUIDE_BINARY_DIR}/WAVE_guide_EN.tex" @ONLY )
   
   add_custom_command (
       OUTPUT ${GUIDE_BINARY_DIR}/WAVE_guide_EN.pdf
       COMMAND ${PDFLATEX_COMPILER} WAVE_guide_EN
       COMMAND ${BIBTEX_COMPILER} WAVE_guide_EN
       COMMAND ${PDFLATEX_COMPILER} WAVE_guide_EN
       COMMAND ${PDFLATEX_COMPILER} WAVE_guide_EN
       DEPENDS ${GUIDE_BINARY_DIR}/WAVE_guide_EN.tex
       WORKING_DIRECTORY ${GUIDE_BINARY_DIR}
   )

   add_custom_target (
       pdf
       DEPENDS ${GUIDE_BINARY_DIR}/WAVE_guide_EN.pdf 
       COMMENT "Creating pdf guide."
   )

else ()

   add_custom_target (
       pdf
       COMMAND echo "ATTENTION: pdflatex/bibtex not found, cannot build user guide" 
   )

endif ()
