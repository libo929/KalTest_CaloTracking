########################################################
# cmake file for building KalFit libraries
# @author Jan Engels, DESY
# @author F. Gaede, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
########################################################


# project name
PROJECT( KalTestCalo )


# project version
SET( KalTest_VERSION_MAJOR 1 )
SET( KalTest_VERSION_MINOR 5 )
SET( KalTest_VERSION_PATCH 4 )

cmake_policy(SET CMP0008 NEW)  

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( ROOT REQUIRED )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )                                                                              
LINK_LIBRARIES( ${ROOT_LIBRARIES} )                                                                                      
ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )  

FIND_PACKAGE( Marlin REQUIRED )
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )                                                                              
LINK_LIBRARIES( ${Marlin_LIBRARIES} )                                                                                      
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )  

FIND_PACKAGE( LCIO REQUIRED ) 
INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )
LINK_LIBRARIES( ${LCIO_LIBRARIES} )
ADD_DEFINITIONS( ${LCIO_DEFINITIONS} )


# macros for generating root dict sources with rootcint
SET( ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
INCLUDE( ${ROOT_DICT_MACROS_FILE} )

### LIBRARY AND TOOLS #######################################################


# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#----- need long long for int64 for now ------
ADD_DEFINITIONS( "-Wno-long-long -fno-strict-aliasing" )


# KalTest library
ADD_SUBDIRECTORY( ./src )


# --------- kaltest examples -------------------------------------------------
# if this option is set to ON examples will be built with default target 'make'
# if set to OFF examples will only be built with target 'make examples'
OPTION( BUILD_EXAMPLES "Set to ON to build examples" OFF )
MESSAGE( STATUS "BUILD_EXAMPLES = ${BUILD_EXAMPLES}" )
ADD_SUBDIRECTORY( ./examples/kaltest )
# ----------------------------------------------------------------------------


#ADD_EXECUTABLE( testcalo src/detector/EXVTXKalTest.cxx )
#TARGET_LINK_LIBRARIES( testcalo KalTestCalo )
#INSTALL( TARGETS testcalo DESTINATION bin )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( KalTestConfig.cmake KalTestConfigVersion.cmake KalTestLibDeps.cmake )

