# Path to tools
export PATH=/prod/tools/common/cmake/3.11.1/bin/:$PATH
export PATH=/prod/tools/common/qt/5.10.1/5.10.1/gcc_64/bin:$PATH
export PATH=/prod/tools/common/qt/5.10.1/Tools/QtCreator/bin:$PATH
export PATH=/prod/tools/rd/python-2.7.10/bin:$PATH

export PATH=/prod/tools/rd/cuda-8.0.61/bin:$PATH

# Library paths
export LD_LIBRARY_PATH=/prod/tools/rd/ffmpeg-2.8.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/prod/tools/rd/libx264/lib:$LD_LIBRARY_PATH


# CMAKE FIND_PACKAGE paths
export CMAKE_PREFIX_PATH=/prod/tools/rd/opencv-3.1.0-noqt:$CMAKE_PREFIX_PATH

# PYTHON
#export PYTHONPATH=/prod/tools/rd/opencv-3.2.0/lib/python2.7/site-packages:$PYTHONPATH

# PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/prod/tools/rd/ffmpeg-2.8.6/lib/pkgconfig
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/prod/tools/rd/libx264/lib/pkgconfig

# Git auto-complete
source /prod/tools/rd/git-2.11.0/share/git-completion.bash

export OpenCV_DIR=/prod/tools/rd/opencv-3.1.0-noqt/share/OpenCV
#export IMAGESLAB=/home/jgagnon/packages/rd_image/1.2.0/libImages.so
#export IMAGESLAB=/prod/tools/rez/packages/int/rd_image/1.1.0/libImages.so


#====================================================================================
source setuphoudini.sh