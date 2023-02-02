#!/bin/bash
 
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 1

#BSUB -J run_beetle_ucv
#BSUB -o run_beetle_ucv.%J.out
#BSUB -e run_beetle_ucv.%J.err 

CURRDIR=$(pwd)
LOGDIRNAME=run_ucv_wind_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_cpu/install/UCV/test_mvgaussian_wind ucv_wind_cpu

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_mvgaussian_wind ucv_wind_gpu

DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/


# openmp case
export OMP_NUM_THREADS=42
export UCV_VTKM_BACKEND=openmp

jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_wind_cpu $DATASETPATH &> ucv_wind_cpu.log

# cuda case
# run two times
export OMP_NUM_THREADS=1

# vtkm only use one gpu
# the time of extracting key will be influenced by the cache time of cpu
# how to reset the gpu memory to avoid it?

export UCV_VTKM_BACKEND=cuda

#double checking what are good proper parameters here
#export UCV_GPU_NUMBLOCK=256
#export UCV_GPU_BLOCKPERTHREAD=128
unset UCV_GPU_NUMBLOCK
unset UCV_GPU_BLOCKPERTHREAD

jsrun -n1 -a1 -c1 -g1 ./ucv_wind_gpu $DATASETPATH  &> ucv_wind_gpu_1.log

jsrun -n1 -a1 -c1 -g1 ./ucv_wind_gpu $DATASETPATH  &> ucv_wind_gpu_2.log


cp *.log $CURRDIR/$LOGDIRNAME

