
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/Matrix.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DTryELEntropy.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DTryELEigenDecomp.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DSampling.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DComputeCases.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DComputeEntropy.hpp"

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

// change this value is number of ensembles changes
constexpr int NumEnsembles = 20;
using SupportedTypesVec = vtkm::List<vtkm::Vec<double, NumEnsembles>>;

void callWorklet(vtkm::cont::Timer &timer, vtkm::cont::DataSet vtkmDataSet, double iso, int numSamples, std::string strategy)
{
  timer.Start();

  vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
  vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
  vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

  // executing the uncertianty thing
  std::cout << "--strategy is " << strategy << "---" << std::endl;

  auto resolveType = [&](const auto &concrete)
  {
    if (strategy == "mvg_multi_worklet")
    {
      vtkm::cont::Invoker invoker;

      // Step 1 compute the eigen decomposition matrix and mean value per cell
      vtkm::cont::ArrayHandle<vtkm::Matrix<vtkm::FloatDefault, 4, 4>> eigenDecomposeMatrix;
      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 4>> meanArray;
      invoker(MVGaussianWithEnsemble2DTryELEigenDecomp{numSamples}, vtkmDataSet.GetCellSet(), concrete, eigenDecomposeMatrix, meanArray);

      std::cout << "step 1 out" << std::endl;
      vtkm::cont::printSummary_ArrayHandle(eigenDecomposeMatrix, std::cout);
      vtkm::cont::printSummary_ArrayHandle(meanArray, std::cout);

      // Step 2 compute sampling array
      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 4>> samplingArray;
      samplingArray.Allocate(numSamples);
      invoker(MVGaussianWithEnsemble2DSampling{}, samplingArray);

      vtkm::cont::printSummary_ArrayHandle(samplingArray, std::cout);

      // Step 3 compute cases for each combination of cell and sample
      vtkm::cont::ArrayHandle<vtkm::UInt8> casesArray;
      // length of casesArray is #cells * # sample
      vtkm::IdComponent numCells = eigenDecomposeMatrix.GetNumberOfValues();
      std::cout << "numCells " << numCells << " numSamples " << numSamples << std::endl;
      casesArray.Allocate(numCells * numSamples);
      invoker(MVGaussianWithEnsemble2DComputeCases{numCells, numSamples, iso}, casesArray, eigenDecomposeMatrix, meanArray, samplingArray);
      
      std::cout << "step 3 out" << std::endl;
      vtkm::cont::printSummary_ArrayHandle(casesArray, std::cout);

      // Step 4 compute corssPorb, entropy, non-zero prob and shrink len of casesArray to #cells
      // there are two ways to do this, one is assign a thread per cell
      // another is assign a thread per sampled data
      crossProbability.Allocate(numCells);
      numNonZeroProb.Allocate(numCells);
      entropy.Allocate(numCells);
      invoker(MVGaussianWithEnsemble2DComputeEntropy{numSamples}, crossProbability, numNonZeroProb, entropy, casesArray);

      // output is cross prob, entropy, num_cross value for each cell position
      vtkm::cont::printSummary_ArrayHandle(crossProbability, std::cout);
    }
    else
    {
      std::cout << "unsupported strategy" << std::endl;
    }
  };

  vtkmDataSet.GetField("ensembles").GetData().CastAndCallForTypes<SupportedTypesVec, VTKM_DEFAULT_STORAGE_LIST>(resolveType);

  timer.Stop();

  // output is ms
  std::cout << "execution time: " << timer.GetElapsedTime() * 1000 << std::endl;

  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << iso;
  std::string isostr = stream.str();

  // check results
  // we use a shallow copy as the data set for
  auto outputDataSet = vtkmDataSet;
  outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);
  outputDataSet.AddCellField("num_nonzero_prob" + isostr, numNonZeroProb);
  outputDataSet.AddCellField("entropy" + isostr, entropy);

  std::string outputFileName = "./test_syntheticdata_el_" + strategy + "_ens" + std::to_string(NumEnsembles) + "_iso" + isostr + ".vtk";
  vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  writeCross.WriteDataSet(outputDataSet);
}

int main(int argc, char *argv[])
{
  vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
      argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
  vtkm::cont::Timer timer{initResult.Device};

  if (argc != 6)
  {
    //./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
    std::cout << "<executable> <SyntheticDataSuffix> <FieldName> <Dim> <iso> <num of sample>" << std::endl;
    exit(0);
  }

  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  std::string dataPathSuffix = std::string(argv[1]);
  std::string fieldName = std::string(argv[2]);

  int dim = std::stoi(argv[3]);
  double isovalue = std::stod(argv[4]);
  int num_samples = std::stoi(argv[5]);

  std::cout << "iso value is: " << isovalue << " num_samples is: " << num_samples << std::endl;

  vtkm::Id xdim = dim;
  vtkm::Id ydim = dim;
  vtkm::Id zdim = 1;

  const vtkm::Id3 dims(xdim, ydim, zdim);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);
  using VecCustomize = vtkm::Vec<double, NumEnsembles>;
  vtkm::cont::ArrayHandle<VecCustomize> dataArraySOA;
  dataArraySOA.Allocate(dim * dim);
  // this is the original data
  // first dim represent number of ensembles
  // second dim represent each version of ensemble for all data points
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> dataArray;

  for (int ensId = 0; ensId < NumEnsembles; ensId++)
  {
    // load each ensemble data
    std::string fileName = dataPathSuffix + "_" + std::to_string(ensId) + ".vtk";
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // aggregate data into one dataset
    // store the ensemble slice by slice
    vtkm::cont::ArrayHandle<vtkm::Float64> fieldDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(fieldName).GetData(), fieldDataArray);

    dataArray.push_back(fieldDataArray);
  }

  // do through the data array and put them into the dataArraySOA
  for (int j = 0; j < dim; j++)
  {
    for (int i = 0; i < dim; i++)
    {
      int index = j * dim + i;

      // each entry has 20 ensembles
      VecCustomize ensembles;
      for (int ensId = 1; ensId <= NumEnsembles; ensId++)
      {
        ensembles[ensId - 1] = dataArray[ensId - 1].ReadPortal().Get(index);
      }

      dataArraySOA.WritePortal().Set(index, ensembles);
    }
  }

  vtkmDataSet.AddPointField("ensembles", dataArraySOA);
  std::cout << "checking input dataset" << std::endl;
  vtkmDataSet.PrintSummary(std::cout);

  callWorklet(timer, vtkmDataSet, isovalue, num_samples, "mvg_multi_worklet");

  return 0;
}