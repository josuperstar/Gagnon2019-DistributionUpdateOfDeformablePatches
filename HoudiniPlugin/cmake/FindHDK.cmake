
#
#    Locate the HDK environment
#
# This module defines
# HDK_FOUND,
# HDK_CXX_COMPILER
# HDK_INCLUDE_DIRS,
# HDK_DEFINITIONS,
# HDK_LIBRARY_DIRS,
# HDK_LIBRARIES,
# HDK_DSO_INSTALL_DIR
#
# For OSX, we have the following as well ...
#
# HDK_FRAMEWORK_DIRS,
# HDK_FRAMEWORKS,
#

SET(HDK_FOUND 1)

#SET(HDK_Version 16.5.439)
SET(HDK_Version 17.0.416)

set(HDK_CXX_COMPILER g++)
set(HDK_INCLUDE_DIRS /prod/software/sidefx/hfs${HDK_Version}/toolkit/include)
#set(HDK_INCLUDE_DIRS /opt/hfs${HDK_Version}/toolkit/include)
# the following are for compiling dso's
set(HDK_DEFINITIONS -DVERSION=\"${HDK_Version}\" -D_GNU_SOURCE -DLINUX -DAMD64 -m64 -fPIC -DSIZEOF_VOID_P=8 -DFBX_ENABLED=1 -DOPENCL_ENABLED=1 -DOPENVDB_ENABLED=1 -DSESI_LITTLE_ENDIAN -DENABLE_THREADS -DUSE_PTHREADS -D_REENTRANT -D_FILE_OFFSET_BITS=64 -c -DGCC4 -DGCC3 -Wno-deprecated -std=c++11 -Wall -W -Wno-parentheses -Wno-sign-compare -Wno-reorder -Wno-uninitialized -Wunused -Wno-unused-parameter -Wno-unused-local-typedefs -DMAKING_DSO -D__UT_DSOVersion__)
set(HDK_LIBRARY_DIRS /usr/X11R6/lib64 /usr/X11R6/lib)
set(HDK_FRAMEWORK_DIRS )
set(HDK_FRAMEWORKS )
set(HDK_LIBRARIES GLU GL X11 Xext Xi dl)
set(HDK_HIH_DIR /home/jo/houdini17.0)


# the following are for compiling dso's
set(HDK_STANDALONE_DEFINITIONS -DVERSION=\"${HDK_Version}\" -D_GNU_SOURCE -DLINUX -DAMD64 -m64 -fPIC -DSIZEOF_VOID_P=8 -DFBX_ENABLED=1 -DOPENCL_ENABLED=1 -DOPENVDB_ENABLED=1 -DSESI_LITTLE_ENDIAN -DENABLE_THREADS -DUSE_PTHREADS -D_REENTRANT -D_FILE_OFFSET_BITS=64 -c -DGCC4 -DGCC3 -Wno-deprecated -std=c++11 -Wall -W -Wno-parentheses -Wno-sign-compare -Wno-reorder -Wno-uninitialized -Wunused -Wno-unused-parameter -Wno-unused-local-typedefs -DMAKING_DSO -D__UT_DSOVersion__)
set(HDK_STANDALONE_LIBRARY_DIRS /prod/software/sidefx/hfs${HDK_Version}/dsolib /usr/X11R6/lib64 /usr/X11R6/lib)
set(HDK_STANDALONE_FRAMEWORK_DIRS )
set(HDK_STANDALONE_FRAMEWORKS )
set(HDK_STANDALONE_LIBRARIES pthread HoudiniUI HoudiniOPZ HoudiniOP3 HoudiniOP2 HoudiniOP1 HoudiniSIM HoudiniGEO HoudiniPRM HoudiniUT GLU GL X11 Xext Xi dl)
