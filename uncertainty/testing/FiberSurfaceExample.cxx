//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
// Example: uniform uncertainty visualization
// Assuming the input is the data that already contains min max array

#include <vtkm/Pair.h>
#include <vtkm/cont/Initialize.h>
#include "../Fiber.h"
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>


int main(int argc, char** argv)
{
  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  vtkm::cont::InitializeResult config = vtkm::cont::Initialize(argc, argv, opts);

  std::string fileName = argv[1];
  std::string outputName = argv[2];
  std::cout << "File Path/File Name" << fileName << std::endl;

  vtkm::io::VTKDataSetReader reader(fileName);
  vtkm::cont::DataSet Data = reader.ReadDataSet();

  vtkm::filter::uncertainty::Fiber filter;
  vtkm::Pair<vtkm::Float64, vtkm::Float64> minAxisValue(0.2, 0.2);
  vtkm::Pair<vtkm::Float64, vtkm::Float64> maxAxisValue(0.3, 0.3);

  filter.SetMaxAxis(maxAxisValue);
  filter.SetMinAxis(minAxisValue);

  filter.SetMinOne("Iron_ensemble_min");
  filter.SetMaxOne("Iron_ensemble_max");
  filter.SetMinTwo("Nickel_ensemble_min");
  filter.SetMaxTwo("Nickel_ensemble_max");

  vtkm::cont::DataSet Output = filter.Execute(Data);
  vtkm::io::VTKDataSetWriter writer(outputName);
  writer.WriteDataSet(Output);

  return 0;
}
