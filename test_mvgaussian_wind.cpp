
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/Initialize.h>

#include "ucvworklet/CreateNewKey.hpp"

#include "ucvworklet/MVGaussianWithEnsemble2DTryLialgEntropy.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DPolyTryLialgEntropy.hpp"

#include <vtkm/filter/clean_grid/CleanGrid.h>
#include <vtkm/filter/geometry_refinement/Triangulate.h>

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

using SupportedTypesVec = vtkm::List<vtkm::Vec<double, 15>>;

void exampleDataSet(int pointNum, std::vector<std::vector<double>> &data)
{

  int version = 15;

  for (int i = 0; i < version; i++)
  {
    std::vector<double> d;
    for (int j = 0; j < pointNum * pointNum; j++)
    {
      d.push_back((j + 0.1) + (i + 1) * j);
    }
    data.push_back(d);
  }
}

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

  // std::string outputFileName = "./red_sea_ucv_iso_" + datatype + "_" + isostr + ".vtk";
  // vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  // writeCross.WriteDataSet(outputDataSet);

  // check results
  // vtkmDataSet.AddCellField("cross_prob", crossProbability);
  std::string outputFileName = "./" + datatype + "_wind_pressure_200_ucv.vtk";
  vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  writeCross.WriteDataSet(outputDataSet);
}

int main(int argc, char *argv[])
{

  vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
      argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
  vtkm::cont::Timer timer{initResult.Device};

  if (argc != 3)
  {
    std::cout << "<executable> <iso> <num of sample>" << std::endl;
    exit(0);
  }

  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  // assuming the ensemble data set is already been extracted out
  // we test results by this dataset
  // https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200
  // to make sure the reuslts are correct

  // load data set, the dim is m*n and for each point there are k ensemble values
  // the data set comes from here https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200

  // vtkm::Id xdim = 240;
  // vtkm::Id ydim = 121;

  vtkm::Id gxdim = 240;
  vtkm::Id gydim = 121;

  // there are some memory issue on cuda if larger than 150*120
  // vtkm::Id xdim = 2;
  // vtkm::Id ydim = 2;
  vtkm::Id xdim = 240;
  vtkm::Id ydim = 121;
  // vtkm::Id zdim = 1;

  // there are cuda errors start with the xdim=40 ydim =40
  // vtkm::Id xdim = 40;
  // vtkm::Id ydim = 40;
  vtkm::Id zdim = 1;

  const vtkm::Id3 dims(xdim, ydim, zdim);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;

  vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);

  std::string windDataDir = "./wind_pressure_200/";
  double isovalue = std::stod(argv[1]);
  int num_samples = std::stoi(argv[2]);
  std::cout << "iso is: " << isovalue << " num_samples is: " << num_samples << std::endl;

  // double isovalue = 1.5;
  //  15 files each contains all data in one file
  //  std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> componentArrays;
  std::vector<std::vector<vtkm::Float64>> dataArray;

  using Vec15 = vtkm::Vec<vtkm::FloatDefault, 15>;

  vtkm::cont::ArrayHandle<Vec15> dataArraySOA;

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> ensemble_min;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> ensemble_max;

  dataArraySOA.Allocate(xdim * ydim);
  ensemble_min.Allocate(xdim * ydim);
  ensemble_max.Allocate(xdim * ydim);

  for (vtkm::IdComponent i = 0; i < 15; i++)
  {
    std::string filename = windDataDir + "Lead_33_" + std::to_string(i) + ".txt";

    // load the data
    std::ifstream fin(filename.c_str());

    std::vector<double> d;
    float element;
    // put data into one vector
    while (fin >> element)
    {
      // put whole picture into the vector
      d.push_back(element);
    }

    // ArrayHandleSOA will group the components in each array and
    // store them into Vec15 automatically.
    // transfer vector into the array handle
    // componentArrays.push_back(vtkm::cont::make_ArrayHandleMove(std::move(d)));
    dataArray.push_back(d);
  }

  // for synthetic dataset
  // exampleDataSet(xdim,dataArray);

  // change aos array to soa array
  // for each points, there are 15 version
  int index = 0;

  // for (vtkm::IdComponent i = 0; i < gxdim * gydim; i++)
  for (vtkm::IdComponent i = 0; i < xdim * ydim; i++)
  {
    Vec15 ensemble;
    vtkm::FloatDefault min = 999999;
    vtkm::FloatDefault max = -999999;
    // if (i == 0 || i == 1 || i == 240 || i == 241)
    {
      // only insert for specific one for testing
      for (vtkm::IdComponent j = 0; j < 15; j++)
      {
        ensemble[j] = dataArray[j][i];
        min = vtkm::Min(min, ensemble[j]);
        max = vtkm::Max(max, ensemble[j]);
      }
      dataArraySOA.WritePortal().Set(index, ensemble);
      ensemble_min.WritePortal().Set(index, min);
      ensemble_max.WritePortal().Set(index, max);

      index++;
    }
  }

  vtkmDataSet.AddPointField("ensembles", dataArraySOA);
  vtkmDataSet.AddPointField("ensembles_min", ensemble_min);
  vtkmDataSet.AddPointField("ensembles_max", ensemble_max);

  // std::cout << "checking input dataset" << std::endl;
  // vtkmDataSet.PrintSummary(std::cout);
  //std::string outputFileName = std::string("wind_pressure_200_MinMax.vtk");
  //vtkm::io::VTKDataSetWriter write(outputFileName);
  //write.WriteDataSet(vtkmDataSet);

  callWorklet(timer, vtkmDataSet, isovalue, num_samples, "stru");
  std::cout << "ok for struc 1" << std::endl;

  callWorklet(timer, vtkmDataSet, isovalue, num_samples, "stru");
  std::cout << "ok for struc 2" << std::endl;

  // test unstructred grid
  // convert the original data to the unstructured grid
  vtkm::filter::clean_grid::CleanGrid clean;
  auto cleanedDataSet = clean.Execute(vtkmDataSet);
  // cleanedDataSet.PrintSummary(std::cout);

  callWorklet(timer, cleanedDataSet, isovalue, num_samples, "unstru");
  std::cout << "ok for unstru" << std::endl;

  vtkm::filter::geometry_refinement::Triangulate triangulate;
  auto tranDataSet = triangulate.Execute(vtkmDataSet);

  // tranDataSet.PrintSummary(std::cout);
  callWorklet(timer, tranDataSet, isovalue, num_samples, "poly");
  std::cout << "ok for poly" << std::endl;

}
