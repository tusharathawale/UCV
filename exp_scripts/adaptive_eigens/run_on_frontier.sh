#!/bin/bash
#SBATCH -A csc143
#SBATCH -J RunAdaptiveEigen
#SBATCH -o %x-%j.out
#SBATCH -t 01:30:00
#SBATCH -p batch
#SBATCH -N 1

DATADIR=/lustre/orion/scratch/zw241/csc331/VisPerfData/uncertainty/
RUNDIR=/lustre/orion/scratch/zw241/csc331/UCVExp
CURRDIR=$(pwd)

mkdir $RUNDIR

cd $RUNDIR

ln -s $CURRDIR/../../install_scripts/frontier_gpu/install/UCV/test_adaptive_eigen test_adaptive_eigen

EIGEN_THRESHOLD=0.01

DIMX=124
DIMY=208
DIMZ=208
NUM_ENS=64
ISO=900
NUM_SAMPLE=500

OMP_THREADS_LIST="64"

for OMP_THREADS in ${OMP_THREADS_LIST}

export OMP_NUM_THREADS=$OMP_THREADS

./test_adaptive_eigen --vtkm-device openmp $DATADIR/beetle_${DIMX}_${DIMY}_${DIMZ}_ens/ens ground_truth $DIMX $DIMY $DIMZ $ISO $NUM_SAMPLE $NUM_ENS beetle_${DIMX}_${DIMY}_${DIMZ}_ens_output $EIGEN_THRESHOLD false &> adaptive_eigen_omp_${OMP_THREADS}.log

done

export OMP_NUM_THREADS=1
# run on hip 
./test_adaptive_eigen --vtkm-device kokkos $DATADIR/beetle_${DIMX}_${DIMY}_${DIMZ}_ens/ens ground_truth $DIMX $DIMY $DIMZ $ISO $NUM_SAMPLE $NUM_ENS beetle_${DIMX}_${DIMY}_${DIMZ}_ens_output $EIGEN_THRESHOLD false &> adaptive_eigen_kokkos.log



