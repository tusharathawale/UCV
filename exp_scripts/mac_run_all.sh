#!/bin/bash
set -x 
set -e

DATADIR=/Users/zhewang/Downloads/UncertaintyCriticalPointDataSets
TESTINGDIR=/Users/zhewang/CodeSpaces/cworkspace/UCV_Critical/UCV/install_scripts/mac/install/UCV

cd $TESTINGDIR

# run on redsea data set
REDSEA_PREFIX=$DATADIR/simplificationRedSeaVTK/simplifiedRedSea
# Run uniform distributino version
./LoadEnsAndProcessByUniform --vtkm-device Any $REDSEA_PREFIX velocityMagnitude 500 500 1 20 &> redsea_uniform.log

# Run MC version, compare each one with the uniform version
./LoadEnsAndProcessByMC --vtkm-device Any $REDSEA_PREFIX velocityMagnitude 500 500 1 20 1000 &> redsea_mc_1000.log

# Run histogram version
./LoadEnsAndProcessByHistogram --vtkm-device Any $REDSEA_PREFIX velocityMagnitude 500 500 1 20 4 &> redsea_bin.log

# Run Epanech version
./LoadEnsAndProcessByEpanech --vtkm-device Any $REDSEA_PREFIX velocityMagnitude 500 500 1 20 &> redsea_epanishkov.log

# Run mv gaussian 
./LoadEnsAndProcessByMVG --vtkm-device Any $REDSEA_PREFIX velocityMagnitude 500 500 1 20 1000 &> redsea_mvg_1000.log

# run on weather data set
WEATHER_DATA=$DATADIR/weather_ens_0121_ens1err_1e-3/ens_0.vtk

./LoadWeatherProcessByUniform --vtkm-device Any $WEATHER_DATA TestField 960 240 1 0.001 &> weather_uniform.log

SAMPLE_NUM=1000
./LoadWeatherProcessByMC --vtkm-device Any $WEATHER_DATA TestField 960 240 1 0.001 $SAMPLE_NUM &> weather_mc_${SAMPLE_NUM}.log

