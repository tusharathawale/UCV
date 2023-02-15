
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/filter/geometry_refinement/Triangulate.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DTryLialgEntropy.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DPolyTryLialgEntropy.hpp"

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using SupportedTypesVec = vtkm::List<vtkm::Vec<double, 20>>;

std::string backend = "openmp";

void initBackend(vtkm::cont::Timer &timer)
{
  // init the vtkh device
  char const *tmp = getenv("UCV_VTKM_BACKEND");

  if (tmp == nullptr)
  {
    std::cout << "no UCV_VTKM_BACKEND env, use openmp" << std::endl;
    backend = "openmp";
  }
  else
  {
    backend = std::string(tmp);
  }

  // if (rank == 0)
  //{
  std::cout << "vtkm backend is:" << backend << std::endl;
  //}

  if (backend == "serial")
  {
    vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
    device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagSerial());
    timer.Reset(vtkm::cont::DeviceAdapterTagSerial());
  }
  else if (backend == "openmp")
  {
    vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
    device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP());
    timer.Reset(vtkm::cont::DeviceAdapterTagOpenMP());
  }
  else if (backend == "cuda")
  {
    vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
    device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagCuda());
    timer.Reset(vtkm::cont::DeviceAdapterTagCuda());
  }
  else
  {
    std::cerr << " unrecognized backend " << backend << std::endl;
  }
  return;
}

void callWorklet(vtkm::cont::Timer &timer, vtkm::cont::DataSet vtkmDataSet, double iso, int numSamples, std::string datatype)
{
  timer.Start();

  vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
  vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
  vtkm::cont::ArrayHandle<vtkm::Float64> entropy;


  if (datatype == "poly")
  {
    // executing the uncertianty thing
    using WorkletType = MVGaussianWithEnsemble2DPolyTryLialgEntropy;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;
    // There are three vertecies
    auto resolveType = [&](const auto &concrete)
    {
      DispatcherType dispatcher(MVGaussianWithEnsemble2DPolyTryLialgEntropy{iso, numSamples});
      dispatcher.Invoke(vtkmDataSet.GetCellSet(), concrete, crossProbability, numNonZeroProb, entropy);
    };

    vtkmDataSet.GetField("ensembles").GetData().CastAndCallForTypes<SupportedTypesVec, VTKM_DEFAULT_STORAGE_LIST>(resolveType);
  }
  else
  {
    // executing the uncertianty thing
    using WorkletType = MVGaussianWithEnsemble2DTryLialgEntropy;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;
    auto resolveType = [&](const auto &concrete)
    {
      DispatcherType dispatcher(MVGaussianWithEnsemble2DTryLialgEntropy{iso, numSamples});
      dispatcher.Invoke(vtkmDataSet.GetCellSet(), concrete, crossProbability, numNonZeroProb, entropy);
    };

    vtkmDataSet.GetField("ensembles").GetData().CastAndCallForTypes<SupportedTypesVec, VTKM_DEFAULT_STORAGE_LIST>(resolveType);
  }

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

  std::string outputFileName = "./red_sea_ucv_iso_" + datatype + "_" + isostr + ".vtk";
  vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  writeCross.WriteDataSet(outputDataSet);
}

static const vtkm::cont::LogLevel CustomLogLevel = vtkm::cont::LogLevel::Perf;

int main(int argc, char *argv[])
{

  if (argc != 3)
  {
    std::cout << "<executable> <iso> <num of sample>" << std::endl;
    exit(0);
  }

  vtkm::cont::SetLogLevelName(CustomLogLevel , "custom");

  vtkm::cont::Initialize(argc, argv);
  vtkm::cont::Timer timer;
  initBackend(timer);
  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  vtkm::Id xdim = 500;
  vtkm::Id ydim = 500;
  vtkm::Id zdim = 1;

  const vtkm::Id3 dims(xdim, ydim, zdim);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;

  vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);

  std::string dataDir = "./red_sea_vtkdata_velocityMagnitude/";
  double isovalue = std::stod(argv[1]);
  int num_samples = std::stoi(argv[2]);
  std::cout << "iso value is: " << isovalue << " num_samples is: " << num_samples << std::endl;

  const int numEnsembles = 20;
  int sliceId = 0;

  using Vec20 = vtkm::Vec<double, numEnsembles>;
  vtkm::cont::ArrayHandle<Vec20> dataArraySOA;
  dataArraySOA.Allocate(xdim * ydim);

  // this is the original data
  // first dim is the ensemble
  // second dim is each slice
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> dataArray;

  for (int ensId = 1; ensId <= numEnsembles; ensId++)
  {
    // load each ensemble data
    std::string fileName = dataDir + "slice_" + std::to_string(sliceId) + "_member_" + std::to_string(ensId) + ".vtk";
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // aggregate data into one dataset
    // store the ensemble slice by slice
    vtkm::cont::ArrayHandle<vtkm::Float64> sliceDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField("velocityMagnitude").GetData(), sliceDataArray);

    dataArray.push_back(sliceDataArray);
  }

  // do through the data array and put them into the dataArraySOA
  for (int j = 0; j < ydim; j++)
  {
    for (int i = 0; i < xdim; i++)
    {
      int index = j * ydim + i;

      // each entry has 20 ensembles
      Vec20 ensembles;
      for (int ensId = 1; ensId <= numEnsembles; ensId++)
      {
        ensembles[ensId - 1] = dataArray[ensId - 1].ReadPortal().Get(index);
      }

      dataArraySOA.WritePortal().Set(index, ensembles);
    }
  }

  vtkmDataSet.AddPointField("ensembles", dataArraySOA);
  std::cout << "checking input dataset" << std::endl;
  vtkmDataSet.PrintSummary(std::cout);

  // std::string outputFileName = "./red_sea_ens_slice_" + std::to_string(sliceId) + ".vtk";
  // vtkm::io::VTKDataSetWriter writeEnsembles(outputFileName);
  // writeEnsembles.WriteDataSet(vtkmDataSet);

  callWorklet(timer, vtkmDataSet, isovalue, num_samples, "stru");
  std::cout << "ok for struc 1" << std::endl;
  
  callWorklet(timer, vtkmDataSet, isovalue, num_samples, "stru");
  std::cout << "ok for struc 2" << std::endl;
  
  // test unstructred grid
  // convert the original data to the unstructured grid
  vtkm::filter::clean_grid::CleanGrid clean;
  auto cleanedDataSet = clean.Execute(vtkmDataSet);
  cleanedDataSet.PrintSummary(std::cout);

  callWorklet(timer, cleanedDataSet, isovalue, num_samples, "unstru");
  std::cout << "ok for unstru" << std::endl;

  vtkm::filter::geometry_refinement::Triangulate triangulate;
  auto tranDataSet = triangulate.Execute(vtkmDataSet);
  tranDataSet.PrintSummary(std::cout);

  callWorklet(timer, tranDataSet, isovalue, num_samples, "poly");
  std::cout << "ok for poly" << std::endl;
  return 0;
}