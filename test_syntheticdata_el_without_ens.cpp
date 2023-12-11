
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

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DTryELEntropy.hpp"

#include "ucvworklet/HelperProbKDE.hpp"
#include "ucvworklet/KDEEntropy.hpp"

#include "ucvworklet/ExtractingMeanStdev.hpp"
#include "ucvworklet/EntropyIndependentGaussian.hpp"

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

void ComputeEntropyWithRuntimeVec(vtkm::cont::DataSet vtkmDataSet,
                                  vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> &runtimeVecArray,
                                  double isovalue, std::string outputFileNameSuffix)
{
  // Processing current ensemble data sets based on uncertianty countour
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdevArray;

  using ComponentType = vtkm::FloatDefault;

  vtkm::cont::Invoker invoke;
  invoke(ExtractingMeanStdevEnsembles{}, runtimeVecArray, meanArray, stdevArray);
  // printSummary_ArrayHandle(meanArray, std::cout);
  // printSummary_ArrayHandle(stdevArray, std::cout);

  vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
  vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
  vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

  invoke(EntropyIndependentGaussian<4, 16>{isovalue}, vtkmDataSet.GetCellSet(), meanArray, stdevArray, crossProbability, numNonZeroProb, entropy);

  auto outputDataSet = vtkmDataSet;

  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << isovalue;
  std::string isostr = stream.str();

  std::string outputFileName = outputFileNameSuffix + isostr + ".vtk";

  outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);
  outputDataSet.AddCellField("num_nonzero_prob" + isostr, numNonZeroProb);
  outputDataSet.AddCellField("entropy" + isostr, entropy);

  vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  writeCross.WriteDataSet(outputDataSet);
}

int main(int argc, char *argv[])
{

  vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
      argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
  vtkm::cont::Timer timer{initResult.Device};

  if (argc != 8)
  {
    //./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
    std::cout << "<executable> <SyntheticDataSuffix> <FieldName> <Dimx> <Dimy> <iso> <num of sampls for mv> <num of ensembles>" << std::endl;
    exit(0);
  }

  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  std::string dataPathSuffix = std::string(argv[1]);
  std::string fieldName = std::string(argv[2]);

  int dimx = std::stoi(argv[3]);
  int dimy = std::stoi(argv[4]);

  double isovalue = std::stod(argv[5]);
  int num_samples = std::stoi(argv[6]);
  int total_num_ensemble = std::stoi(argv[7]);

  std::cout << "iso value is: " << isovalue << " num_samples is: " << num_samples << std::endl;

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
  //for (int ensId = 1; ensId <= total_num_ensemble; ensId++)
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

  std::string outputFileNameSuffix = "./test_2ddata_el_" + fieldName + "_using_all_ens_iso_";
  ComputeEntropyWithRuntimeVec(vtkmDataSet, runtimeVecArray, isovalue, outputFileNameSuffix);
  std::cout << "ok to get entropy for all ensembles" << std::endl;

  // // using ensemble without 0 1 2 ... n-1
  //after loaded, te first element becodes 0
  //for (int noEnsId = 0; noEnsId < total_num_ensemble; noEnsId++)
  for (int noEnsId = 0; noEnsId < total_num_ensemble; noEnsId++)

  {
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> runtimeVecArray(total_num_ensemble - 1);
    runtimeVecArray.Allocate(dimx * dimy);
    auto writePortal = runtimeVecArray.WritePortal();
    for (int j = 0; j < dimy; j++)
    {
      for (int i = 0; i < dimx; i++)
      {
        // for each vertex point, load associated ensemble data
        int pointIndex = j * dimx + i;
        auto vecValue = writePortal.Get(pointIndex);
        // load ensemble data from ens 0 to total-1 except the specific one
        // even though the first ens for red sea is 1 it is 0th element in the loaded array
        int currEnsId = 0;
        int actualEnsId = 0;
        while (actualEnsId < total_num_ensemble)
        {
          if (actualEnsId == noEnsId)
          {
            // skip
            actualEnsId++;
            continue;
          }
          //for debug
          // if (i == 0 && j == 0 && noEnsId==49)
          // {
          //   std::cout << "no ens " << noEnsId << " currEnsId " << currEnsId << " actualEnsId " << actualEnsId << std::endl;
          // }
          vecValue[currEnsId] = dataArray[actualEnsId].ReadPortal().Get(pointIndex);
          actualEnsId++;
          currEnsId++;
        }
      }
    }
    std::string outputFileNameSuffix = "./test_2ddata_el_" + fieldName + "_no_ens_" + std::to_string(noEnsId) + "_iso_";
    ComputeEntropyWithRuntimeVec(vtkmDataSet, runtimeVecArray, isovalue, outputFileNameSuffix);
    std::cout << "ok for no ens " << noEnsId << std::endl;
  }


  return 0;
}

// TODO, maybe than trying to pick out the obvious one and then do similar things again to run same subroutine