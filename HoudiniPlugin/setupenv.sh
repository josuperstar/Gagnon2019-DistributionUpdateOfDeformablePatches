#!/bin/bash
STRING=`pwd`
#source ${ENV_ETC_DIR}/setproject.sh rd_stereo
#settoolset houdini14

echo "[RD SurfaceTextureSynthesis] Setting environment variables to be able to debug and use the local version of the compiled dso."

source /prod/tools/rd/enable
#source /prod/tools/rd/enable-qt5.6.1
export CMAKE_PREFIX_PATH=/prod/tools/rd/opencv-3.1.0-noqt:$CMAKE_PREFIX_PATH
export HOUDINI_DSO_PATH=@/dso_^:@/dso:$STRING/dso/
#export HOUDINI_DSO_PATH=@/dso:$STRING/dso/
echo ${HOUDINI_DSO_PATH}
#export HOUDINI_OTLSCAN_PATH=@/otls_^:@/otls:$STRING/:$HFS/houdini/otls
export HOUDINI_OTLSCAN_PATH=@/otls_^:@/otls:$STRING/DigitalAssets/;$HFS/houdini/otls
echo "Houdini dso and otl have been defined"
cd $STRING


echo "CMAKE_PREFIX_PATH"	+ $CMAKE_PREFIX_PATH
echo "------------------------------------------------"
echo "HOUDINI_DSO_PATH" 	+ $HOUDINI_DSO_PATH
echo "------------------------------------------------"
echo "HOUDINI_OTLSCAN_PATH" 	+ $HOUDINI_OTLSCAN_PATH
echo "------------------------------------------------"

###############################################################################
#				GIT SUBMODULES
#echo "git submodule init and update"

#git submodule init
#git submodule update


#	TEXTURE SYNTHESIS
#cd SurfaceTextureSynthesis
#git submodule init
#git submodule update
#cd $STRING
#echo "git submodule updated"

#	TEXTURE SYNTHESIS
#cd TextureSynthesis
#git submodule init
#git submodule update
#cd $STRING
#echo "git submodule updated"

#	TEXTURE SYNTHESIS
#cd TextureSynthesis/TextureSynthesis
#git submodule init
#git submodule update
#cd $STRING
#echo "git submodule updated"


###############################################################################
#				TEMP DIRECTORIES
# check if temp and render directories exist in HoudiniPlugin/HoudiniScenes/
tempDirectory="./HoudiniScenes/temp"
renderDirectory="./HoudiniScenes/render"
buildDirectory="build"

# bash check if directory exists
if [ -d $renderDirectory ]; then
	echo "Render Directory exists"
else
	echo "Making Render Directory"
	mkdir $renderDirectory
fi

# bash check if directory exists
if [ -d $tempDirectory ]; then
	echo "Temp Directory exists"
else
	echo "Making Temp Directory"
	mkdir $tempDirectory
fi

# bash check if directory exists
if [ -d $buildDirectory ]; then
	echo "Build Directory exists"
else
	echo "Making Build Directory"
	mkdir $buildDirectory
fi
################################################################################


#source /prod/tools/rd/enable
