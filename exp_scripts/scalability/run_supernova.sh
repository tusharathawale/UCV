#!/bin/bash
 
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 43

#BSUB -J run_supernova_log
#BSUB -o run_supernova_log.%J.out
#BSUB -e run_supernova_log.%J.err 

CURRDIR=$(pwd)
LOGDIRNAME=run_supernova_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_mvg_supernova test_mvg_supernova

#soft link for the data
ln -s /gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/supernova_decompose supernova_decompose

export OMP_NUM_THREADS=1
export UCV_VTKM_BACKEND=cuda

#double checking what are good proper parameters here
#export UCV_GPU_NUMBLOCK=256
#export UCV_GPU_BLOCKPERTHREAD=128
unset UCV_GPU_NUMBLOCK
unset UCV_GPU_BLOCKPERTHREAD

PROC_NUM_LIST="4 8 16 32 64 128 256"

for PROC_NUM in ${PROC_NUM_LIST}
do

jsrun -n ${PROC_NUM} -a1 -c1 -g1 ./test_mvg_supernova 0.3 4 1000 false &> test_mvg_supernova_${PROC_NUM}.log

done

cp *.log $CURRDIR/$LOGDIRNAME

