#!/bin/bash
#SBATCH -A csc143
#SBATCH -J FiberUncertainty
#SBATCH -o %x-%j.out
#SBATCH -t 01:30:00
#SBATCH -p batch
#SBATCH -N 1

TestNum=$1
DATADIR=/lustre/orion/scratch/zw241/csc143/VisPerfData/uncertainty/
RUNDIR=/lustre/orion/scratch/zw241/csc143/FiberUncertainty_$TestNum
CURRDIR=$(pwd)

mkdir -p $RUNDIR

cd $RUNDIR

rm TestSuperNova
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/uncertainty/testing/TestSuperNova TestSuperNova
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/uncertainty/testing/TestSuperNovaComparison TestSuperNovaComparison

export OMP_NUM_THREADS=1

# serial

./TestSuperNova --vtkm-device serial $DATADIR/supernova_100_data 64 ClosedForm 5000 &> supernova_100_cf_serial.log

./TestSuperNova --vtkm-device serial $DATADIR/supernova_100_data 64 MonteCarlo 5000 &> supernova_100_mc_serial.log


# kokkos

./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_100_data 64 ClosedForm 5000 &> supernova_100_cf_kokkos.log

./TestSuperNova --vtkm-device kokkos $DATADIR/supernova_100_data 64 MonteCarlo 5000 &> supernova_100_mc_kokkos.log


# openmp
export OMP_NUM_THREADS=64

./TestSuperNova --vtkm-device openmp $DATADIR/supernova_100_data 64 ClosedForm 5000 &> supernova_100_cf_openmp.log

./TestSuperNova --vtkm-device openmp $DATADIR/supernova_100_data 64 MonteCarlo 5000 &> supernova_100_mc_openmp.log


# Testing comparison

./TestSuperNovaComparison --vtkm-device kokkos $DATADIR/supernova_100_data 64 1000 &> supernova_100_comp_samples_1000.log

./TestSuperNovaComparison --vtkm-device kokkos $DATADIR/supernova_100_data 64 2000 &> supernova_100_comp_samples_2000.log

./TestSuperNovaComparison --vtkm-device kokkos $DATADIR/supernova_100_data 64 3000 &> supernova_100_comp_samples_3000.log

./TestSuperNovaComparison --vtkm-device kokkos $DATADIR/supernova_100_data 64 4000 &> supernova_100_comp_samples_4000.log

./TestSuperNovaComparison --vtkm-device kokkos $DATADIR/supernova_100_data 64 5000 &> supernova_100_comp_samples_5000.log

./TestSuperNovaComparison --vtkm-device Serial $DATADIR/supernova_100_data 64 5000 &> supernova_100_comp_samples_5000_serial.log


