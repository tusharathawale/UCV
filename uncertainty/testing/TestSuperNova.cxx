#include <iostream>
#include <random>
#include <vtkm/Pair.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Algorithm.h>
#include "../worklet/ExtractingMinMax.hpp"

#include "../Fiber.h"
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Timer.h>

int main(int argc, char *argv[])
{
  vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
      argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
  std::cout << "initResult.Device: " << initResult.Device.GetName() << std::endl;

  if (argc != 3)
  {
    std::cout << "<executable> <DataFolder> <NumEns>" << std::endl;
    exit(0);
  }

  std::string dataFolder = std::string(argv[1]);
  int NumEns = std::stoi(argv[2]);

  std::string Field1 = "Iron";
  std::string Field2 = "Nickel";

  int dimx, dimy, dimz;

  // load the all ensemble data
  std::vector<vtkm::cont::ArrayHandle<vtkm::FloatDefault>> array1All;
  std::vector<vtkm::cont::ArrayHandle<vtkm::FloatDefault>> array2All;

  for (int i = 0; i < NumEns; i++)
  {
    std::string filePath = dataFolder + "/ens_" + std::to_string(i) + ".vtk";
    vtkm::io::VTKDataSetReader reader(filePath);
    vtkm::cont::DataSet inData = reader.ReadDataSet();
    if (i == 0)
    {
      auto cellSet = inData.GetCellSet();
      vtkm::cont::CellSetStructured<3> structCellSet =
          cellSet.AsCellSet<vtkm::cont::CellSetStructured<3>>();

      vtkm::Id3 pointDims = structCellSet.GetPointDimensions();

      std::cout << "------" << std::endl;
      std::cout << "point dim: " << pointDims[0] << " " << pointDims[1] << " " << pointDims[2] << std::endl;

      // go through all points and set the specific key
      dimx = pointDims[0];
      dimy = pointDims[1];
      dimz = pointDims[2];
    }

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> fieldDataArray1;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(Field1).GetData(), fieldDataArray1);
    array1All.push_back(fieldDataArray1);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> fieldDataArray2;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(Field2).GetData(), fieldDataArray2);
    array2All.push_back(fieldDataArray2);
  }

  std::cout << "cell size " << dimx * dimy * dimz << std::endl;

  // convert the layout, each element is ensemble members
  vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> runtimeVecArray1(NumEns);
  runtimeVecArray1.Allocate(dimx * dimy * dimz);

  vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> runtimeVecArray2(NumEns);
  runtimeVecArray2.Allocate(dimx * dimy * dimz);

  auto writePortal1 = runtimeVecArray1.WritePortal();
  auto writePortal2 = runtimeVecArray2.WritePortal();

  // go through each point, find value, put them into the coresponding ensemble position
  for (int k = 0; k < dimz; k++)
  {
    for (int j = 0; j < dimy; j++)
    {
      for (int i = 0; i < dimx; i++)
      {
        int pointIndex = k * dimx * dimy + j * dimx + i;
        auto vecValue1 = writePortal1.Get(pointIndex);
        auto vecValue2 = writePortal2.Get(pointIndex);

        // load ensemble data from ens 0 to ens with id NumEns-1
        for (int currEnsId = 0; currEnsId < NumEns; currEnsId++)
        {
          // set ensemble value
          vecValue1[currEnsId] = array1All[currEnsId].ReadPortal().Get(pointIndex);
          vecValue2[currEnsId] = array2All[currEnsId].ReadPortal().Get(pointIndex);
        }
      }
    }
  }

  // compute the min and max for each field
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> minField1;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxField1;

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> minField2;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxField2;

  vtkm::cont::Invoker invoke;
  invoke(ExtractingMinMax{}, runtimeVecArray1, minField1, maxField1);
  invoke(ExtractingMinMax{}, runtimeVecArray2, minField2, maxField2);

  // print summary
  // vtkm::cont::printSummary_ArrayHandle(minField1, std::cout);
  // vtkm::cont::printSummary_ArrayHandle(maxField1, std::cout);
  // vtkm::cont::printSummary_ArrayHandle(minField2, std::cout);
  // vtkm::cont::printSummary_ArrayHandle(maxField2, std::cout);

  // user specify the field
  vtkm::filter::uncertainty::Fiber filter;
  vtkm::Pair<vtkm::FloatDefault, vtkm::FloatDefault> minAxisValue(0.2, 0.2);
  vtkm::Pair<vtkm::FloatDefault, vtkm::FloatDefault> maxAxisValue(0.3, 0.3);

  filter.SetMaxAxis(maxAxisValue);
  filter.SetMinAxis(minAxisValue);

  // create the data based on min and max array
  const vtkm::Id3 dims(dimx, dimy, dimz);
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet dataSetForFilter = dataSetBuilder.Create(dims);
  
  //TODO, update spacing based on the original dataset

  dataSetForFilter.AddPointField("ensemble_min_one", minField1);
  dataSetForFilter.AddPointField("ensemble_max_one", maxField1);

  dataSetForFilter.AddPointField("ensemble_min_two", minField2);
  dataSetForFilter.AddPointField("ensemble_max_two", maxField2);

  // call the fiber filter
  filter.SetMinOne("ensemble_min_one");
  filter.SetMaxOne("ensemble_max_one");
  filter.SetMinTwo("ensemble_min_two");
  filter.SetMaxTwo("ensemble_max_two");

  vtkm::cont::Timer timer{initResult.Device};
  std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

  timer.Start();
  vtkm::cont::DataSet output = filter.Execute(dataSetForFilter);
  timer.Stop();
  vtkm::Float64 elapsedTime = timer.GetElapsedTime();
  std::cout << elapsedTime << std::endl;

  //output the dataset
  //vtkm::io::VTKDataSetWriter writer("./out_fiber_supernova_uncertainty_" + std::to_string(dimx) + ".vtk");
  //writer.WriteDataSet(output);

  return 0;
}