#!/bin/bash
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 1

#BSUB -J run_beetle_compare_stage
#BSUB -o run_beetle_compare_stage.%J.out
#BSUB -e run_beetle_compare_stage.%J.err 

CURRDIR=$(pwd)
LOGDIRNAME=run_beetle_compare_stage_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/ucv_reduce_umc ucv_reduce_umc

DATANAME=beetle_124_208_208.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME
FIELD=ground_truth

export OMP_NUM_THREADS=42

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./uvm_point_neighborhood --vtkm-device openmp $DATASETPATH $FIELD mg 4 900 1000 &> ucv_umc_openmp_mg_1000.log

export OMP_NUM_THREADS=1

jsrun -n1 -a1 -c1 -g1 ./uvm_point_neighborhood --vtkm-device cuda $DATASETPATH $FIELD mg 4 900 1000 &> ucv_umc_cuda_mg_1000.log

# kokkos backend

jsrun -n1 -a1 -c1 -g1 ./uvm_point_neighborhood --vtkm-device kokkos $DATASETPATH $FIELD mg 4 900 1000 &> ucv_umc_kokkos_mg_1000.log


# copy things back

cp *.log $CURRDIR/$LOGDIRNAME

# clean the run dir
#rm -r $LOGDIRNAME
