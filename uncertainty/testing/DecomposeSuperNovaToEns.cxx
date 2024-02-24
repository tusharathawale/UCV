// load the original supernova data
// construct the ensemble data from the raw data
// based on spatial downsampling

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

#include "../worklet/ExtractingByNeigoborhoodEnsTwoFields.hpp"

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <filesystem>

using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32,
                                  vtkm::Id>;

int main(int argc, char *argv[])
{

  vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
      argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
  vtkm::cont::Timer timer{initResult.Device};

  if (argc != 4)
  {
    // decompose the 3d data set into the format for a series of ensemble files
    std::cout << "<executable> <FileName> <blocksize> <outputDirName>" << std::endl;
    exit(0);
  }

  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  std::string fileName = std::string(argv[1]);
  std::cout << "fileName: " << fileName << std::endl;

  int blocksize = std::stoi(argv[2]);
  std::string outputDirName = std::string(argv[3]);

  std::string fieldName1 = "Iron";
  std::string fieldName2 = "Nickel";

  // decompose the input into the ensemble data

  vtkm::io::VTKDataSetReader reader(fileName);
  vtkm::cont::DataSet inData = reader.ReadDataSet();

  // check the property of the data
  inData.PrintSummary(std::cout);

  auto field1 = inData.GetField(fieldName1);
  auto field2 = inData.GetField(fieldName2);

  auto cellSet = inData.GetCellSet();

  // Assuming the imput data is the structured data
  bool isStructured = cellSet.IsType<vtkm::cont::CellSetStructured<3>>();
  if (!isStructured)
  {
    std::cout << "the extraction only works for CellSetStructured<3>" << std::endl;
    exit(0);
  }

  vtkm::cont::CellSetStructured<3> structCellSet =
      cellSet.AsCellSet<vtkm::cont::CellSetStructured<3>>();

  vtkm::Id3 pointDims = structCellSet.GetPointDimensions();

  std::cout << "------" << std::endl;
  std::cout << "point dim: " << pointDims[0] << " " << pointDims[1] << " " << pointDims[2] << std::endl;

  vtkm::Id xdim = pointDims[0];
  vtkm::Id ydim = pointDims[1];
  vtkm::Id zdim = pointDims[2];

  vtkm::Id numberBlockx = xdim % blocksize == 0 ? xdim / blocksize : xdim / blocksize + 1;
  vtkm::Id numberBlocky = ydim % blocksize == 0 ? ydim / blocksize : ydim / blocksize + 1;
  vtkm::Id numberBlockz = zdim % blocksize == 0 ? zdim / blocksize : zdim / blocksize + 1;

  // start to extract ensembles
  int ensNum = blocksize * blocksize * blocksize;
  const vtkm::Id3 reducedDims(numberBlockx, numberBlocky, numberBlockz);

  // vtkm introduces some redoundant operations here
  // https://gitlab.kitware.com/vtk/vtk-m/-/merge_requests/3008/pipelines
  auto bounds = inData.GetCoordinateSystem().GetBounds();
  auto reducedOrigin = bounds.MinCorner();
  vtkm::FloatDefault spacex = (bounds.X.Max - bounds.X.Min) / (numberBlockx - 1);
  vtkm::FloatDefault spacey = (bounds.Y.Max - bounds.Y.Min) / (numberBlocky - 1);
  vtkm::FloatDefault spacez = (bounds.Z.Max - bounds.Z.Min) / (numberBlockz - 1);
  vtkm::Vec3f_64 reducedSpaceing(spacex, spacey, spacez);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  // origin is {0,0,0} spacing is {blocksize,blocksize,blocksize} make sure the reduced data
  // are in same shape with original data

  // Check if the directory exists.
  if (std::__fs::filesystem::exists(outputDirName))
  {
    std::__fs::filesystem::remove_all(outputDirName);
  }

  std::__fs::filesystem::create_directory(outputDirName);

  for (vtkm::Id i = 0; i < ensNum; i++)
  {
    // extracting ens id
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> ensValues1;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> ensValues2;
    vtkm::cont::DataSet reducedDataSet = dataSetBuilder.Create(reducedDims, reducedOrigin, reducedSpaceing);
    auto resolveType = [&](const auto &concreteArray1)
    {
      using ArrayType = std::decay_t<decltype(concreteArray1)>;
      auto concreteArray2 = field2.GetData().AsArrayHandle<ArrayType>();
      vtkm::cont::Invoker invoke;
      invoke(ExtractingByNeigoborhoodEnsTwoFields{blocksize, xdim, ydim, zdim, i}, reducedDataSet.GetCellSet(), concreteArray1, ensValues1, concreteArray2, ensValues2);
    };

    field1.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    reducedDataSet.AddPointField(fieldName1, ensValues1);
    reducedDataSet.AddPointField(fieldName2, ensValues2);

    // reducedDataSet.PrintSummary(std::cout);
    // write data into the output dir
    std::string outputFileName = outputDirName + "/ens_" + std::to_string(i) + ".vtk";
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(reducedDataSet);
  }

  return 0;
}

// TODO, maybe than trying to pick out the obvious one and then do similar things again to run same subroutine