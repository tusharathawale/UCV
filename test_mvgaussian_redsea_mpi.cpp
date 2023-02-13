
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

#include <mpi.h>

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

void callWorklet(vtkm::cont::DataSet vtkmDataSet, double iso, int numSamples, std::string datatype)
{

  vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
  vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
  vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

  if (datatype == "poly")
  {

    std::cout << "only support structured case" << std::endl;
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

  /*
  // do not output the data for the parallel case
  // output is ms
  std::cout << "execution time: " << timer.GetElapsedTime() * 1000 << std::endl;

  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << iso;
  std::string isostr = stream.str();

  // check results
  // we use a shallow copy as the data set for
  auto outputDataSet = vtkmDataSet;
  outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);

  std::string outputFileName = "./red_sea_ucv_iso_" + datatype + "_" + isostr + ".vtk";
  vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  writeCross.WriteDataSet(outputDataSet);
  */
}

vtkm::cont::DataSet loadData(int sliceId)
{
  std::string dataDir = "./red_sea_vtkdata_velocityMagnitude/";
  const int numEnsembles = 20;
  vtkm::Id xdim = 500;
  vtkm::Id ydim = 500;
  vtkm::Id zdim = 1;
  using Vec20 = vtkm::Vec<double, numEnsembles>;
  vtkm::cont::ArrayHandle<Vec20> dataArraySOA;
  dataArraySOA.Allocate(xdim * ydim);
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

  const vtkm::Id3 dims(xdim, ydim, zdim);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);
  vtkmDataSet.AddPointField("ensembles", dataArraySOA);
  return vtkmDataSet;
}

int main(int argc, char *argv[])
{
  // compute the number of slices processed by this rank
  int rc = MPI_Init(&argc, &argv);
  int rank;
  int numProcesses;
  if (rc != MPI_SUCCESS)
  {
    printf("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc != 3)
  {
    if (rank == 0)
    {
      std::cout << "<executable> <iso> <num of sample>" << std::endl;
      exit(0);
    }
  }

  double isovalue = std::stod(argv[1]);
  int num_samples = std::stoi(argv[2]);
  if (rank == 0)
  {
    std::cout << "iso value is: " << isovalue << " num_samples is: " << num_samples << std::endl;
  }

  vtkm::cont::Initialize(argc, argv);
  vtkm::cont::Timer timer;
  initBackend(timer);
  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  int totalSlice = 128;
  int actualTotal = 50;
  // there are 50 slices in total
  // for the convenience of performance test
  // use the 128 instead, for the number >=50 slices in total
  // load the existance data again
  if (numProcesses > totalSlice)
  {
    std::cout << "only works when the num of proces < 64" << std::endl;
    exit(0);
  }

  std::vector<vtkm::cont::DataSet> dsList;
  for (int sliceId = 0; sliceId < totalSlice; sliceId++)
  {
    if (sliceId % numProcesses == rank)
    {
      int actualSliceId = sliceId;
      if(actualSliceId>=2*actualTotal){
        actualSliceId = actualSliceId - 2*actualTotal;
      }

      if(actualSliceId>=actualTotal){
        actualSliceId = actualSliceId - actualTotal;
      }

      // current rank load the ith slice
      std::cout << "rank " << rank << " load slice " << actualSliceId << " currid " << sliceId << std::endl;
      vtkm::cont::DataSet ds = loadData(actualSliceId);
      dsList.push_back(ds);
    }
  }

  // do the processing for each member in the list
  // time it
  timer.Start();

  for (int i = 0; i < dsList.size(); i++)
  {
    callWorklet(dsList[i], isovalue, num_samples, "stru");
  }

  //maybe add more operations here
  //such as adding the collective operations based on entropy
  //and compute some results back
  //do some further snalysis here

  MPI_Barrier(MPI_COMM_WORLD);
  timer.Stop();

  if (rank == 0)
  {
    std::cout << "execution time for rank 0: " << timer.GetElapsedTime() * 1000 << std::endl;
  }

  MPI_Finalize();

  return 0;
}