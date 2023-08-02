
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

#include "ucvworklet/HelperComputeDiff.hpp"

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

// There are 20 ensembles for each point
using SupportedTypesVec = vtkm::List<vtkm::Vec<double, 20>, vtkm::Vec<double, 19>>;

void callWorklet(

    vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability,
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb,
    vtkm::cont::ArrayHandle<vtkm::Float64> entropy,
    vtkm::cont::Timer &timer, vtkm::cont::DataSet &vtkmDataSet, double iso, int numSamples, std::string datatype)
{
  timer.Start();

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

  // check results
  // we use a shallow copy as the data set for

  // std::string outputFileName = "./red_sea_ucv_iso_" + datatype + "_" + isostr + "_" + comment +".vtk";
  // vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  // writeCross.WriteDataSet(outputDataSet);
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

  using Vec20 = vtkm::Vec<double, 20>;
  using Vec19 = vtkm::Vec<double, 19>;
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

  // std::string outputFileName = "./red_sea_ens_slice_" + std::to_string(sliceId) + ".vtk";
  // vtkm::io::VTKDataSetWriter writeEnsembles(outputFileName);
  // writeEnsembles.WriteDataSet(vtkmDataSet);
  vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
  vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
  vtkm::cont::ArrayHandle<vtkm::Float64> entropy;
  callWorklet(crossProbability, numNonZeroProb, entropy, timer, vtkmDataSet, isovalue, num_samples, "stru");

  vtkmDataSet.AddCellField("cross_prob_" + std::to_string(isovalue), crossProbability);
  vtkmDataSet.AddCellField("num_nonzero_prob" + std::to_string(isovalue), numNonZeroProb);
  vtkmDataSet.AddCellField("entropy", entropy);

  vtkmDataSet.PrintSummary(std::cout);

  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << isovalue;
  std::string isostr = stream.str();

  std::string outputFileName = "./red_sea_ucv_iso_stru_" + isostr + "_all.vtk";
  vtkm::io::VTKDataSetWriter writeCross(outputFileName);
  writeCross.WriteDataSet(vtkmDataSet);

  // now, we have entropy results for 20 ensemble memebrs
  // we need to execute other 19 times to see what
  // are ensembles without specific slice
  // for (int i = 0; i < numEnsembles; i++)
  for (int nid = 0; nid < numEnsembles; nid++) // try one
  {
    // without ith ensemble
    vtkm::cont::ArrayHandle<Vec19> dataArraySOAInner;
    dataArraySOAInner.Allocate(xdim * ydim);
    for (int j = 0; j < ydim; j++)
    {
      for (int i = 0; i < xdim; i++)
      {
        int index = j * ydim + i;

        // each entry has 20 ensembles
        Vec19 ensemblesInner;
        int count = 0;
        for (int ensId = 0; ensId < numEnsembles; ensId++)
        {
          if (ensId != nid)
          {
            // for each ensemble element, pickout right value and put them into SOA
            ensemblesInner[count] = dataArray[ensId].ReadPortal().Get(index);
            count++;
          }
        }
        dataArraySOAInner.WritePortal().Set(index, ensemblesInner);
      }
    }

    // compute diff between dataArraySOAInner and original dataArraySOA

    vtkm::cont::DataSet vtkmDataSetInner = dataSetBuilder.Create(dims);
    vtkmDataSetInner.AddPointField("ensembles", dataArraySOAInner);
    std::cout << "dataset without ensemble " << nid << std::endl;

    vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

    callWorklet(crossProbability, numNonZeroProb, entropy, timer, vtkmDataSetInner, isovalue, num_samples, "stru");

    vtkmDataSetInner.AddCellField("cross_prob_" + std::to_string(isovalue), crossProbability);
    vtkmDataSetInner.AddCellField("num_nonzero_prob" + std::to_string(isovalue), numNonZeroProb);
    vtkmDataSetInner.AddCellField("entropy", entropy);

    // compute the difference between outputed data set and original one
    // extract the field and compute difference
    // using a simple diff worklet
    vtkm::cont::Invoker invoke;

    // cell data
    vtkm::cont::ArrayHandle<vtkm::Float64> diffArray;
    diffArray.Allocate(vtkmDataSetInner.GetNumberOfCells());

    vtkm::cont::ArrayHandle<vtkm::Float64> field1;
    vtkm::cont::ArrayHandle<vtkm::Float64> field2;

    vtkm::cont::ArrayCopyShallowIfPossible(
        vtkmDataSetInner.GetField("entropy").GetData(), field1);
    vtkm::cont::ArrayCopyShallowIfPossible(
        vtkmDataSet.GetField("entropy").GetData(), field2);

    invoke(HelperComputeDiff{}, field1, field2, diffArray);
    vtkmDataSetInner.AddCellField("entropyDiff", diffArray);

    vtkmDataSetInner.PrintSummary(std::cout);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << isovalue;
    std::string isostr = stream.str();

    std::string outputFileName = "./red_sea_ucv_iso_stru_" + isostr + "_noid" + std::to_string(nid) + ".vtk";
    vtkm::io::VTKDataSetWriter writeCross(outputFileName);
    writeCross.WriteDataSet(vtkmDataSetInner);
  }

  return 0;
}