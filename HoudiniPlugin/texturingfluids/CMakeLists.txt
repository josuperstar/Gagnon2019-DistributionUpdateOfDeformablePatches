# NOTE Change the project name and edit the list of source files to match your project.

# Specify the minimum required version of CMake to build the project
cmake_minimum_required( VERSION 3.0 )

# Set the project name
project( Gagnon2019DistributionUpdateOfDeformablePatches)

##########################
# Initialize the sources #
##########################

# Set the source path
set( SRC_DIR . )

# List the source files
# (except the ones with the Houdini operators)
# NOTE The use of GLOB is not recommended
file( GLOB SOURCES
    ${SRC_DIR}/*.h
    ${SRC_DIR}/*.C
    ${SRC_DIR}/*.cpp
    ${SRC_DIR}/*.hpp
    ${SRC_DIR}/HoudiniPlugins/*
    ${SRC_DIR}/Approaches/*
    ${SRC_DIR}/Strategies/*.[hHcC]*
    ${SRC_DIR}/Core/*.[hHcC]*
    ${SRC_DIR}/Core/Deformations/*
    ${SRC_DIR}/Core/Atlas/*
    ${SRC_DIR}/Core/Gagnon2019/*
    ${SRC_DIR}/PoissonDisk/*
    ${SRC_DIR}/Images/*
    ${SRC_DIR}/Math/*
)

# The output library will have the same name as the project
# and contain all the source files
set( library_name ${PROJECT_NAME} )
add_library( ${library_name} SHARED ${SOURCES} )

# Add include directories
target_include_directories( ${library_name} PRIVATE
    ${SRC_DIR}
    ${Math}
    ${Set}
    ${Images}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# Force c++11
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

##########################
# Register the libraries #
##########################

## Register Houdini (16+) ##

if( DEFINED ENV{HFS} )
    # Display the HFS environment variable
    message( STATUS "HFS directory: " $ENV{HFS} )

    # Locate Houdini's libraries and header files
    list( APPEND CMAKE_PREFIX_PATH "$ENV{HFS}/toolkit/cmake" )
else()
    # If the environment variable HFS is not defined,
    # use -DCMAKE_PREFIX_PATH=/path/to/hfs/toolkit/cmake
    # (or any path to HoudiniConfig.cmake)
    message( STATUS "HFS not set, looking into CMAKE_PREFIX_PATH: " ${CMAKE_PREFIX_PATH} )
    if( NOT CMAKE_PREFIX_PATH )
        message( FATAL_ERROR "Environment variable HFS undefined.
            Please make sure the HFS environment variable is set or
            set CMAKE_PREFIX_PATH to HoudiniConfig.cmake path ({HFS}/toolkit/cmake)." )
    endif()
endif()

# Register an imported library target named 'Houdini'
# (this will use the HoudiniConfig.cmake file)
find_package( Houdini REQUIRED )
message( STATUS "Found Houdini " ${Houdini_VERSION_MAJOR} "." ${Houdini_VERSION_MINOR} "." ${Houdini_VERSION_PATCH} )

# We need version 16 minimum
if( Houdini_VERSION_MAJOR LESS 16 )
    message( FATAL_ERROR "Houdini 16+ required." )
endif()

# Link against the library
target_link_libraries( ${library_name} Houdini )

# Include ${CMAKE_CURRENT_BINARY_DIR} for the generated headers
target_include_directories( ${library_name} PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

## Register OpenCV ##

# In Windows, the OPENCV_DIR environment variable should have been set
# (see https://docs.opencv.org/master/d3/d52/tutorial_windows_install.html)
# and can be used to find OpenCVConfig.cmake
if( DEFINED ENV{OPENCV_DIR} )
    message( STATUS "OpenCV directory: " $ENV{OPENCV_DIR} )
    list( APPEND CMAKE_PREFIX_PATH "$ENV{OPENCV_DIR}/lib" )
endif()
find_package( OpenCV REQUIRED PATHS "/Volumes/Projects/Production/packages/tools/opencv/4.1.2.2/bin")

# Link against the libraries
target_link_libraries( ${library_name} ${OpenCV_LIBS} )

# Add the necessary include directories
target_include_directories( ${library_name} PRIVATE ${OpenCV_INCLUDE_DIRS} )

#################
# Houdini steps #
#################

# Configure several common target properties
# We leave the default (Houdini's default dso directory)
#houdini_configure_target( ${library_name} )
# Otherwise, we can set the output directory (INSTDIR) manually
#houdini_configure_target( ${library_name} INSTDIR ${CMAKE_CURRENT_LIST_DIR}/dso )
houdini_configure_target(${library_name} INSTDIR ${CMAKE_CURRENT_LIST_DIR}/../../setup/rez/dso)
