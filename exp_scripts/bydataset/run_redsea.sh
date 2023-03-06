#!/bin/bash
 
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 1

#BSUB -J run_redsea
#BSUB -o run_redsea.%J.out
#BSUB -e run_redsea.%J.err 

# This scrip need to be updated 

CURRDIR=$(pwd)
LOGDIRNAME=run_redsea_log

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_mvgaussian_redsea test_mvgaussian_redsea

ln -s /gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/red_sea_vtkdata_velocityMagnitude red_sea_vtkdata_velocityMagnitude

# openmp case
export OMP_NUM_THREADS=42

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_redsea --vtkm-device=openmp 0.1 1000 &> openmp_1000.log
 
jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_redsea --vtkm-device=openmp 0.1 2000 &> openmp_2000.log

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./test_mvgaussian_redsea --vtkm-device=openmp 0.1 4000 &> openmp_4000.log


export OMP_NUM_THREADS=1

#serial adaptor
jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=serial 0.1 1000 &> serial_1000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=serial 0.1 2000 &> serial_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=serial 0.1 4000 &> serial_4000.log


# cuda case
# run two times
jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=cuda 0.1 1000 &> cuda_1000_1.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=cuda 0.1 1000 &> cuda_1000_2.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=cuda 0.1 2000 &> cuda_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=cuda 0.1 4000 &> cuda_4000.log

#kokkos cuda case
export OMP_NUM_THREADS=1
jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 1000 &> kokkos_1000_1.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 1000 &> kokkos_1000_2.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 2000 &> kokkos_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 4000 &> kokkos_4000.log

cp *.log $CURRDIR/$LOGDIRNAME

