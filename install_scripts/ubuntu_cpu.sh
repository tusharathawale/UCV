#!/bin/bash
set -e

build_jobs=6
mkdir -p ubuntu_cpu
cd ubuntu_cpu

HERE=`pwd`
source $HERE/../settings.sh

SOFTWARE_SRC_DIR="$HERE/src"
SOFTWARE_BUILD_DIR="$HERE/build"
SOFTWARE_INSTALL_DIR="$HERE/install"

mkdir -p $SOFTWARE_SRC_DIR
mkdir -p $SOFTWARE_BUILD_DIR
mkdir -p $SOFTWARE_INSTALL_DIR



echo "====> Installing EasyLinalg"
EASY_LINALG_SRC_DIR="$SOFTWARE_SRC_DIR/EasyLinalg"
EASY_LINALG_INSTALL_DIR="$HERE/../../ucvworklet/linalg/EasyLinalg/"

rm -rf $EASY_LINALG_SRC_DIR
cd $SOFTWARE_SRC_DIR
git clone $EASY_LINALG_REPO

# move include dir to correct place

# clean old dir if it exist
if [ -d $EASY_LINALG_INSTALL_DIR ]; then
    rm -rf $EASY_LINALG_INSTALL_DIR
fi

mkdir -p $EASY_LINALG_INSTALL_DIR

# move files to new dir
cp EasyLinalg/StaticMemTemplate/include/* $EASY_LINALG_INSTALL_DIR
# clean source files
rm -rf $EASY_LINALG_SRC_DIR

echo "====> Installing EasyLinalg, ok"

echo "====> build UCV"
# the only have build dir without the install dir
# the install dir is same with the build dir
UCV_SRC_DIR=$HERE/../../
# use the install dir as the build dir
UCV_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/UCV"

#if [ -d $UCV_INSTALL_DIR ]; then
    echo "====> skip, $UCV_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
#else
    
    #-DCMAKE_BUILD_TYPE=Release \
    # using Debug mode for enabling the assert option in testing
    cmake -B ${UCV_INSTALL_DIR} -S ${UCV_SRC_DIR} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_PARAVIEW_PLUGIN=ON \
    -DParaView_DIR=/home/zw1/Documents/cworkspace/build_paraview_serial/lib/cmake/paraview-5.12 \
    -DVTKm_DIR=/home/zw1/Documents/cworkspace/build_paraview_serial/lib/cmake/paraview-5.12/vtk/vtkm 

    cd $HERE

    # build and install
    echo "**** Building UCV"
    cmake --build ${UCV_INSTALL_DIR} -j${build_jobs}
#fi

# not sure why the libvtkmdiympi.so is not included during the build process
echo "try to add library path by executing:"
echo "export LD_LIBRARY_PATH=${VTKM_INSTALL_DIR}/lib:\${LD_LIBRARY_PATH}"
