#!/bin/bash
#SBATCH -A csc143
#SBATCH -J FiberUncertainty
#SBATCH -o %x-%j.out
#SBATCH -t 01:30:00
#SBATCH -p batch
#SBATCH -N 1

TestNum=$1
DATADIR=/lustre/orion/scratch/zw241/csc143/VisPerfData/uncertainty/
RUNDIR=/lustre/orion/scratch/zw241/csc143/FiberUncertainty_RedSea
CURRDIR=$(pwd)

mkdir -p $RUNDIR

cd $RUNDIR

rm TestRedSea
rm TestRedSeaComparison
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/uncertainty/testing/TestRedSea TestRedSea
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/uncertainty/testing/TestRedSeaComparison TestRedSeaComparison

export OMP_NUM_THREADS=1

# serial
./TestRedSea --vtkm-device serial $DATADIR/redsea &> redsea_serial.log

# kokkos

./TestRedSea --vtkm-device kokkos $DATADIR/redsea &> redsea_kokkos.log

# openmp
export OMP_NUM_THREADS=64

./TestRedSea --vtkm-device openmp $DATADIR/redsea &> redsea_openmp.log

# Testing comparison

./TestRedSeaComparison --vtkm-device kokkos $DATADIR/redsea 1000 &> redsea_comp_samples_1000.log

./TestRedSeaComparison --vtkm-device kokkos $DATADIR/redsea 2000 &> redsea_comp_samples_2000.log

./TestRedSeaComparison --vtkm-device kokkos $DATADIR/redsea 3000 &> redsea_comp_samples_3000.log

./TestRedSeaComparison --vtkm-device kokkos $DATADIR/redsea 4000 &> redsea_comp_samples_4000.log

./TestRedSeaComparison --vtkm-device kokkos $DATADIR/redsea 5000 &> redsea_comp_samples_5000.log

./TestRedSeaComparison --vtkm-device serial $DATADIR/redsea 5000 &> redsea_comp_samples_serial_5000.log
