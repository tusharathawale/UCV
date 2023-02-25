#!/bin/bash
set -e

module load rocm/5.3.0
module load cmake
module load git-lfs
#module load eigen/3.3.9

build_jobs=4
mkdir -p crusher_gpu
cd crusher_gpu

HERE=`pwd`
source $HERE/../settings.sh

SOFTWARE_SRC_DIR="$HERE/src"
SOFTWARE_BUILD_DIR="$HERE/build"
SOFTWARE_INSTALL_DIR="$HERE/install"

mkdir -p $SOFTWARE_SRC_DIR
mkdir -p $SOFTWARE_BUILD_DIR
mkdir -p $SOFTWARE_INSTALL_DIR

echo "====> Install Kokkos"
kokkos_src_dir=${SOFTWARE_SRC_DIR}/kokkos
kokkos_build_dir=${SOFTWARE_BUILD_DIR}/kokkos
kokkos_install_dir=${SOFTWARE_INSTALL_DIR}/kokkos

if false; then
rm -rf ${kokkos_src_dir}
git clone -b master https://github.com/kokkos/kokkos.git ${kokkos_src_dir}
rm -rf ${kokkos_build_dir}
cmake -S ${kokkos_src_dir} -B ${kokkos_build_dir} \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON\
  -DKokkos_ARCH_VEGA90A=ON \
  -DCMAKE_CXX_COMPILER=hipcc \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=OFF \
  -DCMAKE_INSTALL_PREFIX=${kokkos_install_dir} \
  -DCMAKE_CXX_FLAGS="--amdgpu-target=gfx90a"
cmake --build ${kokkos_build_dir} -j10
cmake --install ${kokkos_build_dir}
fi

echo "====> Installing vtk-m"
VTKM_SRC_DIR="${SOFTWARE_SRC_DIR}/vtk-m"
VTKM_BUILD_DIR="${SOFTWARE_BUILD_DIR}/vtk-m"
VTKM_INSTALL_DIR="${SOFTWARE_INSTALL_DIR}/vtk-m"

# check the install dir
if [ -d $VTKM_INSTALL_DIR ]; then
    echo "====> skip, $VTKM_INSTALL_DIR already exists," \
             "please remove it if you want to reinstall it"
else
    echo $VTKM_SRC_DIR
    echo $VTKM_BUILD_DIR
    echo $VTKM_INSTALL_DIR
    # check vktm source dir
    if [ ! -d $VTKM_SRC_DIR ]; then
    # clone the source
    cd $SOFTWARE_SRC_DIR
    git clone $VTKM_REPO
    cd $VTKM_SRC_DIR
    git checkout master
    fi
    
    cd $HERE

    # build and install
    echo "**** Building vtk-m"

    # TODO, the gpu version can be different here
    # we only use the cpu version here
    # there are still some issues to run gpu and cpu backend
    # by the same binary? the gpu is dorced to be used anyway?

    cmake -B ${VTKM_BUILD_DIR} -S ${VTKM_SRC_DIR} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=OFF \
    -DVTKm_USE_DEFAULT_TYPES_FOR_ASCENT=ON \
    -DVTKm_USE_DOUBLE_PRECISION=ON \
    -DVTKm_NO_DEPRECATED_VIRTUAL=ON \
    -DVTKm_USE_64BIT_IDS=OFF \
    -DCMAKE_INSTALL_PREFIX=${VTKM_INSTALL_DIR} \
    -DVTKm_ENABLE_MPI=ON \
    -DVTKm_ENABLE_OPENMP=ON \
    -DVTKm_ENABLE_LOGGING=ON \
    -DVTKm_ENABLE_RENDERING=ON \
    -DVTKm_ENABLE_KOKKOS=ON \
    -DVTKm_ENABLE_TESTING=OFF \
    -DCMAKE_HIP_ARCHITECTURES=gfx90a \
    -DCMAKE_PREFIX_PATH=${kokkos_install_dir} \
    -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_C_COMPILER=hipcc

    cmake --build ${VTKM_BUILD_DIR} -j${build_jobs}

    echo "**** Installing vtk-m"
    cmake --install ${VTKM_BUILD_DIR}
fi

echo "====> Installing vtk-m, ok"

echo "vtkm install dir ${VTKM_INSTALL_DIR}"

echo "====> build UCV"
# the only have build dir without the install dir
# the install dir is same with the build dir
UCV_SRC_DIR=$HERE/../../
# use the install dir as the build dir
UCV_INSTALL_DIR="$SOFTWARE_INSTALL_DIR/UCV"

#if [ -d $UCV_INSTALL_DIR ]; then
#    echo "====> skip, $UCV_INSTALL_DIR already exists," \
#             "please remove it if you want to reinstall it"
#else

    cmake -B ${UCV_INSTALL_DIR} -S ${UCV_SRC_DIR} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=OFF \
    -DUSE_HIP=ON \
    -DCMAKE_PREFIX_PATH=${kokkos_install_dir} \
    -DVTKm_DIR=${VTKM_INSTALL_DIR}/lib/cmake/vtkm-2.0 \
    -DCMAKE_HIP_ARCHITECTURES=gfx90a \
    -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_C_COMPILER=hipcc
    
    cd $HERE

    # build and install
    echo "**** Building UCV"
    cmake --build ${UCV_INSTALL_DIR} -j${build_jobs}
#fi

# not sure why the libvtkmdiympi.so is not included during the build process
echo "try to add library path by executing:"
echo "export LD_LIBRARY_PATH=${VTKM_INSTALL_DIR}/lib:\${LD_LIBRARY_PATH}"
