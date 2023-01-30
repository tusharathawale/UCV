#!/bin/bash
 
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 1

#BSUB -J run_beetle_ucv
#BSUB -o run_beetle_ucv.%J.out
#BSUB -e run_beetle_ucv.%J.err 

CURRDIR=$(pwd)
LOGDIRNAME=run_beetle_ucv_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_cpu/install/UCV/ucv_reduce_umc ucv_umc_cpu

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/ucv_reduce_umc ucv_umc_gpu

DATANAME=beetle_496_832_832.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME

# set openmp thread

export OMP_NUM_THREADS=1

jsrun -n1 -a1 -c1 -g0 -bpacked:1 ./ucv_umc_cpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_cpu_serial_uni.log

jsrun -n1 -a1 -c1 -g0 -bpacked:1 ./ucv_umc_cpu $DATASETPATH ground_truth ig 4 900 &> ucv_umc_cpu_serial_ig.log

jsrun -n1 -a1 -c1 -g0 -bpacked:1 ./ucv_umc_cpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_cpu_serial_mg.log


export OMP_NUM_THREADS=42

export UCV_VTKM_BACKEND=openmp

jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_umc_cpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_cpu_uni.log

jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_umc_cpu $DATASETPATH ground_truth ig 4 900 &> ucv_umc_cpu_ig.log

jsrun -n1 -a1 -c42 -g0 -bpacked:42 ./ucv_umc_cpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_cpu_mg.log


export OMP_NUM_THREADS=1

# vtkm only use one gpu
# the time of extracting key will be influenced by the cache time of cpu
# how to reset the gpu memory to avoid it?

export UCV_VTKM_BACKEND=cuda

jsrun -n1 -a1 -c1 -g1 ./ucv_umc_gpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_gpu_uni_1.log

# executing it another time, it seems there are some extra work or cache for the first tun
# the uniform distribution
jsrun -n1 -a1 -c1 -g1 ./ucv_umc_gpu $DATASETPATH ground_truth uni 4 900 &> ucv_umc_gpu_uni_2.log


# the indepedent gaussin distribution
jsrun -n1 -a1 -c1 -g1 ./ucv_umc_gpu $DATASETPATH ground_truth ig 4 900 &> ucv_umc_gpu_ig.log


# the multivariant gaussin distribution
jsrun -n1 -a1 -c1 -g1 ./ucv_umc_gpu $DATASETPATH ground_truth mg 4 900 &> ucv_umc_gpu_mg.log

# add the results for the mg 

# copy things back

cp *.log $CURRDIR/$LOGDIRNAME

# clean the run dir
#cd ..

#rm -r $LOGDIRNAME
