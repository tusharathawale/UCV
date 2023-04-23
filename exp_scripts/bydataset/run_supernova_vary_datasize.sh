#!/bin/bash
#BSUB -P csc143
#BSUB -W 01:59
#BSUB -nnodes 1

#BSUB -J run_supernova_vary_datasize
#BSUB -o run_supernova_vary_datasize.%J.out
#BSUB -e run_supernova_vary_datasize.%J.err 

CURRDIR=$(pwd)
DATANAME_LIST="supernova_resample_100.vtk supernova_resample_200.vtk supernova_resample_400.vtk supernova_resample_800.vtk"

for DATANAME in ${DATANAME_LIST}
do

LOGDIRNAME=run_supernova_vary_datasize_log_${DATANAME}

mkdir -p $CURRDIR/$LOGDIRNAME

cd $MEMBERWORK/csc143

rm -rf $LOGDIRNAME
mkdir $LOGDIRNAME

cd $LOGDIRNAME

ln -s $CURRDIR/../../install_scripts/summit_gpu/install/UCV/ucv_reduce_umc ucv_reduce_umc

DATASETPATH=/gpfs/alpine/proj-shared/csc143/zhewang/datasets/uncertainty/$DATANAME
FIELD=Iron
ISO=0.3

export OMP_NUM_THREADS=42

jsrun -n1 -a1 -c42 -g1 -bpacked:42 ./ucv_reduce_umc --vtkm-device openmp $DATASETPATH $FIELD mg 4 $ISO 1000 &> ucv_umc_openmp_mg_1000.log

export OMP_NUM_THREADS=1

jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device cuda $DATASETPATH $FIELD mg 4 $ISO 1000 &> ucv_umc_cuda_mg_1000.log

# kokkos backend

jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device kokkos $DATASETPATH $FIELD mg 4 $ISO 1000 &> ucv_umc_kokkos_mg_1000.log

# serial backend
jsrun -n1 -a1 -c1 -g1 ./ucv_reduce_umc --vtkm-device serial $DATASETPATH $FIELD mg 4 $ISO 1000 &> ucv_umc_serial_mg_1000.log

cp *.log $CURRDIR/$LOGDIRNAME

# clean the run dir
#rm -r $LOGDIRNAME

done