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

# DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/300_300_one_peak_gaussian_sigma_0.02/RawdataPointScalar
DATASETPATH_200=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/200_200_move_xy_one_pattern/RawdataPointScalar
DATASETPATH_400=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/400_400_move_xy_one_pattern/RawdataPointScalar
DATASETPATH_800=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/800_800_move_xy_one_pattern/RawdataPointScalar

FIELD=TestField

export OMP_NUM_THREADS=42

#For openmp backend
jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_syntheticdata_el --vtkm-device openmp $DATASETPATH_200 $FIELD 200 0.8 1000 &> one_worklet_el_openmp_200_1000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_el_multiple_worklet --vtkm-device openmp $DATASETPATH_200 $FIELD 200 0.8 1000 &> multi_worklet_openmp_200_1000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_syntheticdata_el --vtkm-device openmp $DATASETPATH_200 $FIELD 200 0.8 2000 &> one_worklet_el_openmp_200_2000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_el_multiple_worklet --vtkm-device openmp $DATASETPATH_200 $FIELD 200 0.8 2000 &> multi_worklet_openmp_200_2000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_syntheticdata_el --vtkm-device openmp $DATASETPATH_400 $FIELD 400 0.8 1000 &> one_worklet_el_openmp_400_1000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_el_multiple_worklet --vtkm-device openmp $DATASETPATH_400 $FIELD 400 0.8 1000 &> multi_worklet_openmp_400_1000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_syntheticdata_el --vtkm-device openmp $DATASETPATH_400 $FIELD 400 0.8 2000 &> one_worklet_el_openmp_400_2000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_el_multiple_worklet --vtkm-device openmp $DATASETPATH_400 $FIELD 400 0.8 2000 &> multi_worklet_openmp_400_2000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_syntheticdata_el --vtkm-device openmp $DATASETPATH_800 $FIELD 800 0.8 1000 &> one_worklet_el_openmp_800_1000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_el_multiple_worklet --vtkm-device openmp $DATASETPATH_800 $FIELD 800 0.8 1000 &> multi_worklet_openmp_800_1000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_syntheticdata_el --vtkm-device openmp $DATASETPATH_800 $FIELD 800 0.8 2000 &> one_worklet_el_openmp_800_2000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_el_multiple_worklet --vtkm-device openmp $DATASETPATH_800 $FIELD 800 0.8 2000 &> multi_worklet_openmp_800_2000.log

# testing cuda backend

export OMP_NUM_THREADS=1

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH_200 $FIELD 200 0.8 1000 &> one_worklet_el_cuda_200_1000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH_200 $FIELD 200 0.8 1000 &> multi_worklet_cuda_200_1000.log

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH_200 $FIELD 200 0.8 2000 &> one_worklet_el_cuda_200_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH_200 $FIELD 200 0.8 2000 &> multi_worklet_cuda_200_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH_400 $FIELD 400 0.8 1000 &> one_worklet_el_cuda_400_1000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH_400 $FIELD 400 0.8 1000 &> multi_worklet_cuda_400_1000.log

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH_400 $FIELD 400 0.8 2000 &> one_worklet_el_cuda_400_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH_400 $FIELD 400 0.8 2000 &> multi_worklet_cuda_400_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH_800 $FIELD 800 0.8 1000 &> one_worklet_el_cuda_800_1000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH_800 $FIELD 800 0.8 1000 &> multi_worklet_cuda_800_1000.log

jsrun -n1 -a1 -c1 -g1 ./test_syntheticdata_el --vtkm-device cuda $DATASETPATH_800 $FIELD 800 0.8 2000 &> one_worklet_el_cuda_800_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_el_multiple_worklet --vtkm-device cuda $DATASETPATH_800 $FIELD 800 0.8 2000 &> multi_worklet_cuda_800_2000.log

# it seems that the multi-worklet way is not efficient as the one worklet approach

cp *.log $CURRDIR/$LOGDIRNAME
