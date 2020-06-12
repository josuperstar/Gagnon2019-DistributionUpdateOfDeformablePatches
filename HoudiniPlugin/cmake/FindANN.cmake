# This module finds libann.
#
# It sets the following variables:
#  ANN_FOUND              - Set to false, or undefined, if lemon isn't found.
#  ANN_INCLUDE_DIR        - include directory.
#  ANN_LIBRARIES          - library files

SET(ANN_INCLUDE_DIR /prod/tools/rd/ann1.1.2/include)
SET(ANN_LIBRARIES /prod/tools/rd/ann1.1.2/lib/libANN.so)

#FIND_PATH(ANN_INCLUDE_DIR ${ANN_INCLUDE}/ANN.h PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include ENV C_INCLUDE_PATH ENV CPLUS_INCLUDE_PATH )
#FIND_PATH(ANN_INCLUDE_DIR /automount/sun-01/home/jgagnon/mokko/dltdev/HoudiniPlugin/VolumeGenerator/ANN/include/ANN )
FIND_LIBRARY(ANN_LIBRARIES ann PATHS ENV LD_LIBRARY_PATH ENV LIBRARY_PATH)

GET_FILENAME_COMPONENT(ANN_LIBRARY_PATH $ANN_LIBRARIES} PATH)
SET( ANN_LIBRARY_DIR ${LEMON_LIBRARY_PATH} CACHE PATH "Path to ANN library.")

# handle the QUIETLY and REQUIRED arguments and set ANN_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ANN DEFAULT_MSG ANN_LIBRARIES ANN_INCLUDE_DIR)
#FIND_PACKAGE_HANDLE_STANDARD_ARGS(ANN DEFAULT_MSG ANN_INCLUDE_DIR)

#MARK_AS_ADVANCED( ANN_INCLUDE_DIR ANN_LIBRARIES ANN_LIBRARY_DIR )