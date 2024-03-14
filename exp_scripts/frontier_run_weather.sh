#!/bin/bash
#SBATCH -A csc143
#SBATCH -J CriticalPoint
#SBATCH -o %x-%j.out
#SBATCH -t 01:30:00
#SBATCH -p batch
#SBATCH -N 1

TestNum=$1
DATADIR=/lustre/orion/scratch/zw241/csc143/VisPerfData/uncertainty/weather_data/weather_ens_0121_ens1err_1e-3/ens_0.vtk
RUNDIR=/lustre/orion/scratch/zw241/csc143/CriticalPoint_$TestNum
CURRDIR=$(pwd)

mkdir -p $RUNDIR

cd $RUNDIR

rm LoadWeatherProcessByMC
rm LoadWeatherProcessByUniform
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/LoadWeatherProcessByMC LoadWeatherProcessByMC
ln -s $CURRDIR/../install_scripts/frontier_gpu/install/UCV/LoadWeatherProcessByUniform LoadWeatherProcessByUniform

export OMP_NUM_THREADS=1

# The Uniform version

./LoadWeatherProcessByUniform --vtkm-device kokkos $DATADIR TestField 960 240 1 0.001 &> weather_uniform.log

# Run MC version, compare each one with the uniform version
SAMPLE_NUM_LIST="100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000"

for SAMPLE_NUM in ${SAMPLE_NUM_LIST}
do

./LoadWeatherProcessByMC --vtkm-device kokkos $DATADIR TestField 960 240 1 0.001 $SAMPLE_NUM &> weather_mc_${SAMPLE_NUM}.log

python3 $CURRDIR/../compute_diff.py MinProb_MC_Weather960_240_${SAMPLE_NUM}.vtk MinProb_Uniform_Weather960_240.vtk MinProb &> diff_${SAMPLE_NUM}.log 

done

