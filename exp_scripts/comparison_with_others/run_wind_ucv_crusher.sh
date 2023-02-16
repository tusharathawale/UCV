#!/bin/bash
 
#SBATCH -A CSC331
#SBATCH -J test
#SBATCH -o test.out
#SBATCH -t 01:00:00
#SBATCH -N 1

# This scrip need to be updated 

CURRDIR=$(pwd)
LOGDIRNAME=run_ucv_wind_log_crusher_$1

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

#ln -s $CURRDIR/../../install_scripts/summit_cpu/install/UCV/test_mvgaussian_wind ucv_wind_cpu

ln -s $CURRDIR/../../install_scripts/crusher_gpu/install/UCV/test_mvgaussian_wind ucv_wind_gpu

DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/wind_pressure_200 

cp -r $DATASETPATH .


# running example on original run in debug node
# jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./run ../datasets/txt_files/

# openmp case
export OMP_NUM_THREADS=42
export UCV_VTKM_BACKEND=openmp

#jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_wind_cpu &> ucv_wind_cpu.log

# cuda case
# run two times
export OMP_NUM_THREADS=1

# vtkm only use one gpu
# the time of extracting key will be influenced by the cache time of cpu
# how to reset the gpu memory to avoid it?

export UCV_VTKM_BACKEND=hip

#double checking what are good proper parameters here
#export UCV_GPU_NUMBLOCK=256
#export UCV_GPU_BLOCKPERTHREAD=128
unset UCV_GPU_NUMBLOCK
unset UCV_GPU_BLOCKPERTHREAD


srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_wind_gpu 0.3 1000 &> ucv_wind_gpu_1.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_wind_gpu 0.3 1000 &> ucv_wind_gpu_2.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_wind_gpu 0.3 2000 &> ucv_wind_gpu_3.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_wind_gpu 0.3 4000 &> ucv_wind_gpu_4.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_wind_gpu 0.3 8000 &> ucv_wind_gpu_5.log

cp *.log $CURRDIR/$LOGDIRNAME

