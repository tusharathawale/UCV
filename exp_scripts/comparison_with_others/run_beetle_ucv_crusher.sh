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
LOGDIRNAME=run_beetle_ucv_log_crusher_$1

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

#ln -s $CURRDIR/../../install_scripts/summit_cpu/install/UCV/ucv_reduce_umc ucv_umc_cpu

ln -s $CURRDIR/../../install_scripts/crusher_gpu/install/UCV/ucv_reduce_umc ucv_umc_gpu

DATANAME=beetle_496_832_832.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME

# set openmp thread

export OMP_NUM_THREADS=1

#jsrun -n1 -a1 -c1 -g0 -bpacked:1 ./ucv_umc_cpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_cpu_serial_uni.log

#jsrun -n1 -a1 -c1 -g0 -bpacked:1 ./ucv_umc_cpu $DATASETPATH ground_truth ig 4 900 &> ucv_umc_cpu_serial_ig.log

#this one is time consuming, add it when it is necessary
#jsrun -n1 -a1 -c1 -g0 -bpacked:1 ./ucv_umc_cpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_cpu_serial_mg.log


export OMP_NUM_THREADS=42

export UCV_VTKM_BACKEND=openmp

#jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_umc_cpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_cpu_uni.log

#jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_umc_cpu $DATASETPATH ground_truth ig 4 900 &> ucv_umc_cpu_ig.log

#jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_umc_cpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_cpu_mg.log


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

#./ucv_umc_gpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_gpu_uni_1.log

#./ucv_umc_gpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_gpu_uni_2.log

#./ucv_umc_gpu $DATASETPATH ground_truth ig 4 900 &> ucv_umc_gpu_ig.log

#./ucv_umc_gpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_gpu_mg.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_gpu_uni_1.log

# executing it another time, it seems there are some extra work or cache for the first tun
# the uniform distribution
srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_gpu_uni_2.log


# the indepedent gaussin distribution
srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu $DATASETPATH ground_truth ig 4 900 &> ucv_umc_gpu_ig.log


# the multivariant gaussin distribution
srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_gpu_mg.log

# add the results for the mg 

# copy things back

cp *.log $CURRDIR/$LOGDIRNAME

# clean the run dir
#cd ..

#rm -r $LOGDIRNAME
