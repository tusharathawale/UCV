
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

// using 1 ensmeble, compute entropy
// using 2 ensemble, compute entropy
// ...
// using 15 ensembles
// compute how info changes using compute_entropy_change.py

constexpr int NumEnsembles = 7;
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
    if (strategy == "ig")
    {
      // extracting mean and stdev
      vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
      vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdevArray;

      using DispatcherType = vtkm::worklet::DispatcherMapField<ExtractingMeanStdevEnsembles>;
      DispatcherType dispatcher;
      dispatcher.Invoke(concrete, meanArray, stdevArray);

      using WorkletType = EntropyIndependentGaussian<4, 16>;
      using DispatcherEntropyIG = vtkm::worklet::DispatcherMapTopology<WorkletType>;

      DispatcherEntropyIG dispatcherEntropyIG(EntropyIndependentGaussian<4, 16>{iso});
      dispatcherEntropyIG.Invoke(vtkmDataSet.GetCellSet(), meanArray, stdevArray, crossProbability, numNonZeroProb, entropy);
    }
    else if (strategy == "mvg")
    {
      using WorkletType = MVGaussianWithEnsemble2DTryELEntropy;
      using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;
      DispatcherType dispatcher(MVGaussianWithEnsemble2DTryELEntropy{iso, numSamples});
      dispatcher.Invoke(vtkmDataSet.GetCellSet(), concrete, crossProbability, numNonZeroProb, entropy);
    }
    else if (strategy == "kde")
    {

      // compute kde distribution
      using DispatcherProbKDE = vtkm::worklet::DispatcherMapField<HelperProbKDE<NumEnsembles>>;

      DispatcherProbKDE dispatcherProbKDE(HelperProbKDE<NumEnsembles>{iso});
      vtkm::cont::ArrayHandle<vtkm::FloatDefault> postiveProb;

      dispatcherProbKDE.Invoke(concrete, postiveProb);

      // compute entropy things
      // go through cell by points
      using DispatcherEntropyKDE = vtkm::worklet::DispatcherMapTopology<KDEEntropy<4, 16>>;

      DispatcherEntropyKDE dispatcherEntropyKDE(KDEEntropy<4, 16>{iso});
      dispatcherEntropyKDE.Invoke(vtkmDataSet.GetCellSet(), postiveProb, crossProbability, numNonZeroProb, entropy);
    }
    else if (strategy == "mvkde")
    {
      // TODO
      // need to do the similar things with kde
      // adjust the positive probabilty based on density value
      // ? How to decide the positiv value, 1d case has the elf function
      // should we also do the sampling for this case?
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

  std::string outputFileName = "./test_syntheticdata_el_sequence_" + strategy + "_ens" + std::to_string(NumEnsembles) + "_iso" + isostr + ".vtk";
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
  //redsea data start from 1
  //for (int ensId = 1; ensId <= total_num_ensemble; ensId++)
  for (int ensId = 0; ensId < total_num_ensemble; ensId++)
  {
    std::string fileName = dataPathSuffix + "_" + std::to_string(ensId) + ".vtk";
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::Float64> fieldDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(fieldName).GetData(), fieldDataArray);
    dataArray.push_back(fieldDataArray);
  }

  // start from using 1 ensemble to all ensembles
  // for each ensemble value, extract corresponding values from dataArray
  for (int usedEnsembles = 1; usedEnsembles <= total_num_ensemble; usedEnsembles++)
  // for (int usedEnsembles = 1; usedEnsembles <= 3; usedEnsembles++)
  {
    // init runtime vec
    // in each point, there are usedEnsembles values
    // there are dimx*dimy points in total for 2d data sets
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> runtimeVecArray(usedEnsembles);
    runtimeVecArray.Allocate(dimx * dimy);
    auto writePortal = runtimeVecArray.WritePortal();
    for (int j = 0; j < dimy; j++)
    {
      for (int i = 0; i < dimx; i++)
      {
        int pointIndex = j * dimx + i;
        auto vecValue = writePortal.Get(pointIndex);

        // load ensemble data from ens 0 to ens with id usedEnsembles-1
        for (int currEndId = 0; currEndId < usedEnsembles; currEndId++)
        {
          // set ensemble value
          vecValue[currEndId] = dataArray[currEndId].ReadPortal().Get(pointIndex);
        }
      }
    }

    // checking results
    std::cout << "ok for loading data with " << usedEnsembles << " ensembles, start processing it" << std::endl;
    // printSummary_ArrayHandle(runtimeVecArray, std::cout);

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

    outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);
    outputDataSet.AddCellField("num_nonzero_prob" + isostr, numNonZeroProb);
    outputDataSet.AddCellField("entropy" + isostr, entropy);

    std::string outputFileName = "./test_2ddata_el_" + fieldName + "_using_ens_" + std::to_string(usedEnsembles) + "_iso" + isostr + ".vtk";
    vtkm::io::VTKDataSetWriter writeCross(outputFileName);
    writeCross.WriteDataSet(outputDataSet);
  }

  return 0;
}

//TODO, maybe than trying to pick out the obvious one and then do similar things again to run same subroutine