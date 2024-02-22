//
//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <iostream>
#include <random>
#include <vtkm/Pair.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/filter/uncertainty/Fiber.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

namespace
{
template <typename T>
vtkm::cont::DataSet MakeLogValuesTestDataSet()
{
  const vtkm::Id3 dims(25, 25, 25);

  vtkm::Id numPoints = dims[0] * dims[1] * dims[2];
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet dataSet = dataSetBuilder.Create(dims);

  std::vector<T> ensemble_max_one;
  std::vector<T> ensemble_min_one;
  std::vector<T> ensemble_max_two;
  std::vector<T> ensemble_min_two;
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<double> minValue(-20, -1);
  std::uniform_real_distribution<double> maxValue(0, 20);
  for (vtkm::Id i = 0; i < numPoints; ++i)
  {
    double randomMinValueOne = minValue(gen);
    double randomMaxValueOne = maxValue(gen);
    double randomMinValueTwo = minValue(gen);
    double randomMaxValueTwo = maxValue(gen);

    ensemble_max_one.push_back(static_cast<T>(randomMaxValueOne));
    ensemble_min_one.push_back(static_cast<T>(randomMinValueOne));
    ensemble_min_two.push_back(static_cast<T>(randomMinValueTwo));
    ensemble_max_two.push_back(static_cast<T>(randomMaxValueTwo));
  }

  dataSet.AddPointField("ensemble_max_one", ensemble_max_one);
  dataSet.AddPointField("ensemble_min_one", ensemble_min_one);
  dataSet.AddPointField("ensemble_max_two", ensemble_max_two);
  dataSet.AddPointField("ensemble_min_two", ensemble_min_two);
  return dataSet;
}

void TestUncertaintyGeneral()
{

  vtkm::cont::DataSet input = MakeLogValuesTestDataSet<vtkm::FloatDefault>();
  vtkm::filter::uncertainty::Fiber filter;

  vtkm::Pair<vtkm::Float64, vtkm::Float64> minAxisValue(0.2, 0.2);
  vtkm::Pair<vtkm::Float64, vtkm::Float64> maxAxisValue(0.205, 0.205);

  filter.SetMaxAxis(maxAxisValue);
  filter.SetMinAxis(minAxisValue);
  filter.SetMinOne("ensemble_min_one");
  filter.SetMaxOne("ensemble_max_one");
  filter.SetMinTwo("ensemble_max_two");
  filter.SetMaxTwo("ensemble_min_two");

  vtkm::cont::DataSet output = filter.Execute(input);

  vtkm::io::VTKDataSetWriter writer("/Users/n5j/Desktop/out_fiber_uncertainty.vtk");
  writer.WriteDataSet(output);
}

void TestFiber()
{
  TestUncertaintyGeneral();
}
} // anonymous namespace

int Fiber(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(TestFiber, argc, argv);
}
