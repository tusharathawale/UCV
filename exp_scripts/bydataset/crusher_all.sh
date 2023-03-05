#!/bin/bash

#SBATCH -A CSC331
#SBATCH -J test
#SBATCH -o test.out
#SBATCH -t 01:59:00
#SBATCH -N 1

ml cmake rocm/5.3.0

set -e
set -x

CURRDIR=$(pwd)
LOGDIRNAME=crusher_all_$1

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

#ln -s $CURRDIR/../../install_scripts/summit_cpu/install/UCV/ucv_reduce_umc ucv_umc_cpu
ln -s $CURRDIR/../../install_scripts/crusher_gpu/install/UCV/ucv_reduce_umc ucv_umc_gpu
ln -s $CURRDIR/../../install_scripts/crusher_gpu/install/UCV/uvm_point_neighborhood uvm_point_neighborhood

DATANAME=beetle_496_832_832.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME

export OMP_NUM_THREADS=1

# using the reduced by key
srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu --vtkm-device kokkos $DATASETPATH ground_truth uni 4 900 &> beetle_uni_1.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu --vtkm-device kokkos $DATASETPATH ground_truth uni 4 900 &> beetle_uni_2.log
# the indepedent gaussin distributio
srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu --vtkm-device kokkos $DATASETPATH ground_truth ig 4 900 &> beetle_ig.log

# the multivariant gaussin distribution, this does not work on crusher
# srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./ucv_umc_gpu --vtkm-device kokkos $DATASETPATH ground_truth mg 4 900 &> ucv_umc_gpu_mg.log

# using the fixed block size
srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./uvm_point_neighborhood --vtkm-device kokkos $DATASETPATH ground_truth uni 4 900 &> beetle_point_neighborhood_uni_1.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./uvm_point_neighborhood --vtkm-device kokkos $DATASETPATH ground_truth uni 4 900 &> beetle_point_neighborhood_uni_2.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./uvm_point_neighborhood --vtkm-device kokkos $DATASETPATH ground_truth ig 4 900 &> beetle_point_neighborhood_ig.log


# using the smaller data size, check beetle mg

DATANAME=beetle_124_208_208.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME

jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device kokkos $DATASETPATH $FIELD mg 4 900 1000 &> beetle_small_mg_1000.log


# checking the supernova
DATANAME=supernova_visit_400_400_400.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME
FIELD=Iron
ISO=0.3

jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device kokkos $DATASETPATH $FIELD uni 4 $ISO 1000 &> super_uni_1.log

jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device kokkos $DATASETPATH $FIELD uni 4 $ISO 1000 &> super_uni_2.log

jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device kokkos $DATASETPATH $FIELD ig 4 $ISO 1000 &> super_ig.log


# using the smaller data size, check supernova mg

DATANAME=supernova_visit_152_152_152.vtk
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME

jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device kokkos $DATASETPATH Iron mg 4 0.3 1000 &> supernova_small_mg_1000.log

# wind data set
ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_mvgaussian_wind test_mvgaussian_wind
DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/wind_pressure_200 

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 1000 &> wind_mg_1000_1.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 1000 &> wind_mg_1000_2.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 2000 &> wind_mg_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 4000 &> wind_mg_4000.log


# red sea data set
ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/test_mvgaussian_redsea test_mvgaussian_redsea
ln -s /gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/red_sea_vtkdata_velocityMagnitude red_sea_vtkdata_velocityMagnitude

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 1000 &> redsea_mg_1000_1.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 1000 &> redsea_mg_1000_2.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 2000 &> redsea_mg_2000.log

jsrun -n1 -a1 -c1 -g1 ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 4000 &> redsea_mg_4000.log


# copy things back


cp *.log $CURRDIR/$LOGDIRNAME

# clean the run dir
#cd ..

#rm -r $LOGDIRNAME