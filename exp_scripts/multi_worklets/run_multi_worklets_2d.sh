#!/bin/bash
#BSUB -P csc143
#BSUB -W 01:00
#BSUB -nnodes 1

#BSUB -J run_multi_worklets_2d
#BSUB -o run_multi_worklets_2d.%J.out
#BSUB -e run_multi_worklets_2d.%J.err 

CURRDIR=$(pwd)
LOGDIRNAME=run_multi_worklets_2d_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_syntheticdata_el test_syntheticdata_el
ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_mvgaussian_el_multiple_worklet test_mvgaussian_el_multiple_worklet

DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/300_300_one_peak_gaussian_sigma_0.02/RawdataPointScalar
FIELD=TestField

export OMP_NUM_THREADS=42

#there are issues if we set g as 0 even if for the openmp backend
jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_syntheticdata_el --vtkm-device openmp $DATASETPATH $FIELD 300 0.8 1000 &> test_syntheticdata_el_openmp.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_el_multiple_worklet --vtkm-device openmp $DATASETPATH $FIELD 300 0.8 1000 &> multi_worklet_openmp.log

export OMP_NUM_THREADS=1

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH $FIELD 300 0.8 1000 &> test_syntheticdata_el_openmp_cuda_1.log

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH $FIELD 300 0.8 1000 &> test_syntheticdata_el_openmp_cuda_2.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH $FIELD 300 0.8 1000 &> multi_worklet_cuda_1.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH $FIELD 300 0.8 1000 &> multi_worklet_cuda_2.log

# it seems that the multi-worklet way is not efficient as the one worklet approach

cp *.log $CURRDIR/$LOGDIRNAME
