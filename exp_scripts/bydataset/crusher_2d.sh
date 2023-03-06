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
LOGDIRNAME=crusher_2d_$1

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/crusher_gpu/install/UCV/test_mvgaussian_wind test_mvgaussian_wind
# wind data set
ln -s /gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/wind_pressure_200 wind_pressure_200

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 1000 &> wind_mg_1000_1.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 1000 &> wind_mg_1000_2.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 2000 &> wind_mg_2000.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_wind --vtkm-device=kokkos 0.3 4000 &> wind_mg_4000.log


# red sea data set
ln -s $CURRDIR/../../install_scripts/crusher_gpu/install/UCV/test_mvgaussian_redsea test_mvgaussian_redsea
ln -s /gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/red_sea_vtkdata_velocityMagnitude red_sea_vtkdata_velocityMagnitude

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 1000 &> redsea_mg_1000_1.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 1000 &> redsea_mg_1000_2.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 2000 &> redsea_mg_2000.log

srun -A CSC331 -N 1 -n 1 --ntasks-per-node=1 --gpus-per-node=1 --gpu-bind=closest ./test_mvgaussian_redsea --vtkm-device=kokkos 0.1 4000 &> redsea_mg_4000.log


# copy things back


cp *.log $CURRDIR/$LOGDIRNAME

# clean the run dir
#cd ..

#rm -r $LOGDIRNAME