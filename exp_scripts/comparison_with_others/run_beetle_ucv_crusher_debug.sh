#!/bin/bash

#SBATCH -A CSC331
#SBATCH -J test
#SBATCH -o test.out
#SBATCH -t 00:15:00
#SBATCH -N 1

ml cmake rocm/5.3.0

set -e
set -x

CURRDIR=$(pwd)
LOGDIRNAME=run_beetle_ucv_log_crusher_debug_$1

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/crusher_gpu/install/UCV/ucv_reduce_umc ucv_umc_gpu

DATANAME=beetle_8_16_16.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME

export OMP_NUM_THREADS=1
export UCV_VTKM_BACKEND=hip

# the multivariant gaussin distribution
srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_gpu_mg.log

