#!/bin/bash
 
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 9

#BSUB -J run_redsea_log
#BSUB -o run_redsea_log.%J.out
#BSUB -e run_redsea_log.%J.err 

CURRDIR=$(pwd)
LOGDIRNAME=run_redsea_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_mvgaussian_redsea_mpi test_mvgaussian_redsea_mpi

#soft link for the data
ln -s /gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/red_sea_vtkdata_velocityMagnitude red_sea_vtkdata_velocityMagnitude

export OMP_NUM_THREADS=1
export UCV_VTKM_BACKEND=cuda

#double checking what are good proper parameters here
#export UCV_GPU_NUMBLOCK=256
#export UCV_GPU_BLOCKPERTHREAD=128
unset UCV_GPU_NUMBLOCK
unset UCV_GPU_BLOCKPERTHREAD

PROC_NUM_LIST="2 4 8 16 32 50"

for PROC_NUM in ${PROC_NUM_LIST}
do

jsrun -n ${PROC_NUM} -a1 -c1 -g1 ./test_mvgaussian_redsea_mpi 0.1 1000 &> ucv_redsea_proc_${PROC_NUM}.log

done

cp *.log $CURRDIR/$LOGDIRNAME

