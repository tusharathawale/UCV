#!/bin/bash
#SBATCH -A csc143
#SBATCH -J CriticalPoint
#SBATCH -o %x-%j.out
#SBATCH -t 01:30:00
#SBATCH -p batch
#SBATCH -N 1

TestNum=$1
DATADIR=/lustre/orion/scratch/zw241/csc143/VisPerfData/uncertainty/red_sea_vtkdata_velocityMagnitude_slice10/slice_10_member
RUNDIR=/lustre/orion/scratch/zw241/csc143/CriticalPoint_$TestNum
CURRDIR=$(pwd)

mkdir -p $RUNDIR

cd $RUNDIR

rm LoadEnsAndProcessByUniform
rm LoadEnsAndProcessByMC
rm LoadEnsAndProcessByHistogram
rm LoadEnsAndProcessByEpanech

ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/LoadEnsAndProcessByUniform LoadEnsAndProcessByUniform
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/LoadEnsAndProcessByMC LoadEnsAndProcessByMC
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/LoadEnsAndProcessByHistogram LoadEnsAndProcessByHistogram
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/LoadEnsAndProcessByEpanech LoadEnsAndProcessByEpanech

export OMP_NUM_THREADS=1

# The Uniform version
./LoadEnsAndProcessByUniform --vtkm-device kokkos $DATADIR velocityMagnitude 500 500 1 20 &> redsea_uniform.log

# Run MC version, compare each one with the uniform version
./LoadEnsAndProcessByMC --vtkm-device kokkos $DATADIR velocityMagnitude 500 500 1 20 1000 &> redsea_mc_1000.log

# Run histogram version
./LoadEnsAndProcessByHistogram --vtkm-device kokkos $DATADIR velocityMagnitude 500 500 1 20 4 &> redsea_bin.log

# Run Epanech version
./LoadEnsAndProcessByEpanech --vtkm-device kokkos $DATADIR velocityMagnitude 500 500 1 20 &> redsea_epanishkov.log