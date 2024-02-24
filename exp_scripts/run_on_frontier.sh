#!/bin/bash
#SBATCH -A csc143
#SBATCH -J FiberUncertainty
#SBATCH -o %x-%j.out
#SBATCH -t 01:30:00
#SBATCH -p batch
#SBATCH -N 1

DATADIR=/lustre/orion/scratch/zw241/csc143/VisPerfData/uncertainty/
RUNDIR=/lustre/orion/scratch/zw241/csc143/FiberUncertainty
CURRDIR=$(pwd)

mkdir -p $RUNDIR

cd $RUNDIR

rm TestSuperNova
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/uncertainty/testing/TestSuperNova TestSuperNova

export OMP_NUM_THREADS=1

# serial
./TestSuperNova --vtkm-device serial $DATADIR/supernova_25_data 64 &> supernova_25_serial.log

./TestSuperNova --vtkm-device serial $DATADIR/supernova_50_data 64 &> supernova_50_serial.log

./TestSuperNova --vtkm-device serial $DATADIR/supernova_100_data 64 &> supernova_100_serial.log

# kokkos
./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_25_data 64 &> supernova_25_kokkos_1.log

./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_25_data 64 &> supernova_25_kokkos_2.log

./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_50_data 64 &> supernova_50_kokkos_1.log

./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_50_data 64 &> supernova_50_kokkos_2.log

./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_100_data 64 &> supernova_100_kokkos_1.log

./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_100_data 64 &> supernova_100_kokkos_2.log

# openmp
export OMP_NUM_THREADS=64

./TestSuperNova --vtkm-device openmp $DATADIR/supernova_25_data 64 &> supernova_25_openmp.log

./TestSuperNova --vtkm-device openmp $DATADIR/supernova_50_data 64 &> supernova_50_openmp.log

./TestSuperNova --vtkm-device openmp $DATADIR/supernova_100_data 64 &> supernova_100_openmp.log

