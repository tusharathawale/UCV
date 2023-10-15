
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
#include "ucvworklet/MVGaussianWithEnsemble2DTryELEntropy.hpp"

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using SupportedTypesVec = vtkm::List<vtkm::Vec<double, 20>>;

void callWorklet(vtkm::cont::Timer &timer, vtkm::cont::DataSet vtkmDataSet,double iso, int numSamples, std::string datatype)
{
  timer.Start();

  vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
  vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
  vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

  if (datatype == "poly")
  {
    std::cout << "poly is not suppoted yet" << std::endl;
  }
  else
  {
    // executing the uncertianty thing
    std::cout << "--Test using easy liag library---" << std::endl;
    using WorkletType = MVGaussianWithEnsemble2DTryELEntropy;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;
    auto resolveType = [&](const auto &concrete)
    {
      DispatcherType dispatcher(MVGaussianWithEnsemble2DTryELEntropy{iso, numSamples});
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

  std::string outputFileName = "./test_syntheticdata_el_" + datatype + "_" + isostr + ".vtk";
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

  const int numEnsembles = 20;

  vtkm::Id xdim = dim;
  vtkm::Id ydim = dim;
  vtkm::Id zdim = 1;

  const vtkm::Id3 dims(xdim, ydim, zdim);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);
  using Vec20 = vtkm::Vec<double, numEnsembles>;
  vtkm::cont::ArrayHandle<Vec20> dataArraySOA;
  dataArraySOA.Allocate(dim * dim);
  // this is the original data
  // first dim represent number of ensembles
  // second dim represent each version of ensemble for all data points
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> dataArray;

  for (int ensId = 0; ensId < numEnsembles; ensId++)
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

  /*

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
  */
  return 0;
}