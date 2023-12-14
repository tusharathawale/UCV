
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandleRuntimeVec.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

int main(int argc, char *argv[])
{

  vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
      argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
  vtkm::cont::Timer timer{initResult.Device};

  if (argc != 6)
  {
    //./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
    std::cout << "<executable> <SyntheticDataSuffix> <FieldName> <Dimx> <Dimy> <num of ensembles>" << std::endl;
    exit(0);
  }

  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  std::string dataPathSuffix = std::string(argv[1]);
  std::string fieldName = std::string(argv[2]);

  int dimx = std::stoi(argv[3]);
  int dimy = std::stoi(argv[4]);

  int total_num_ensemble = std::stoi(argv[5]);

  vtkm::Id xdim = dimx;
  vtkm::Id ydim = dimy;
  vtkm::Id zdim = 1;

  const vtkm::Id3 dims(xdim, ydim, zdim);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);

  // load all data values
  vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> allEnsemblesArray(total_num_ensemble);
  allEnsemblesArray.Allocate(dimx * dimy);

  // for each ensemble(version) of the data
  // store results into the allEnsemblesArray
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> dataArray;
  // redsea data start from 1
  for (int ensId = 0; ensId < total_num_ensemble; ensId++)
  // for (int ensId = 1; ensId <= total_num_ensemble; ensId++)
  {
    std::string fileName = dataPathSuffix + "_" + std::to_string(ensId) + ".vtk";
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::Float64> fieldDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(fieldName).GetData(), fieldDataArray);
    dataArray.push_back(fieldDataArray);
  }

  std::cout << "ok to load the data at the first step" << std::endl;

  // using all ensembles
  vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> runtimeVecArray(total_num_ensemble);
  runtimeVecArray.Allocate(dimx * dimy);
  auto writePortal = runtimeVecArray.WritePortal();
  for (int j = 0; j < dimy; j++)
  {
    for (int i = 0; i < dimx; i++)
    {
      int pointIndex = j * dimx + i;
      auto vecValue = writePortal.Get(pointIndex);
      // load ensemble data from ens 0 to ens with id usedEnsembles-1
      for (int currEndId = 0; currEndId < total_num_ensemble; currEndId++)
      {
        // set ensemble value
        vecValue[currEndId] = dataArray[currEndId].ReadPortal().Get(pointIndex);
      }
    }
  }

  std::string outputFileNameSuffix = "./ensemble_data_" + fieldName + "_ens" + std::to_string(total_num_ensemble) + ".vtk";
  vtkmDataSet.AddPointField("ensemble", runtimeVecArray);
  vtkm::io::VTKDataSetWriter writeCross(outputFileNameSuffix);
  writeCross.WriteDataSet(vtkmDataSet);

  return 0;
}

// TODO, maybe than trying to pick out the obvious one and then do similar things again to run same subroutine