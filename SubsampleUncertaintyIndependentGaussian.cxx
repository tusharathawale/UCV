//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "SubsampleUncertaintyIndependentGaussian.h"

#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayHandleUniformPointCoordinates.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ErrorBadType.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/Keys.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/ExtractingMeanStdev.hpp"

namespace
{

VTKM_CONT void ComputeMeanStdevForField(
    const vtkm::filter::uncertainty::SubsampleUncertaintyIndependentGaussian* self,
    vtkm::cont::DataSet& data,
    const vtkm::cont::Field& field,
    const vtkm::worklet::Keys<vtkm::Id>& keys)
{
  if (field.IsPointField())
  {
    vtkm::cont::Invoker invoke;
    vtkm::cont::UnknownArrayHandle inArray = field.GetData();
    vtkm::cont::UnknownArrayHandle meanArray = inArray.NewInstanceFloatBasic();
    vtkm::cont::UnknownArrayHandle stdevArray = inArray.NewInstanceFloatBasic();
    // These allocations are not necessary in the most recent version of VTK-m.
    meanArray.Allocate(keys.GetInputRange());
    stdevArray.Allocate(keys.GetInputRange());
    auto resolveType = [&](const auto& concrete) {
      auto meanConcrete = meanArray.ExtractArrayFromComponents<vtkm::FloatDefault>();
      auto stdevConcrete = stdevArray.ExtractArrayFromComponents<vtkm::FloatDefault>();
      invoke(ExtractingMeanStdev{}, keys, concrete, meanConcrete, stdevConcrete);
    };
    inArray.CastAndCallWithExtractedArray(resolveType);
    data.AddPointField(field.GetName() + self->GetMeanSuffix(), meanArray);
    data.AddPointField(field.GetName() + self->GetStdevSuffix(), stdevArray);
  }
  else if (field.IsWholeDataSetField())
  {
    // No change.
    data.AddField(field);
  }
  else
  {
    // Cell fields not supported. They get dropped.
  }
}

} // anonymous namespace

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

vtkm::cont::DataSet SubsampleUncertaintyIndependentGaussian::DoExecute(const vtkm::cont::DataSet& input)
{
  if (!input.GetCellSet().IsType<vtkm::cont::CellSetStructured<3>>())
  {
    throw vtkm::cont::ErrorBadType("Extraction only works for CellSetStructured<3>.");
  }
  vtkm::cont::CellSetStructured<3> cellSet;
  input.GetCellSet().AsCellSet(cellSet);

  vtkm::Id3 numPoints = cellSet.GetPointDimensions();

  vtkm::Id3 numBlocks = (numPoints + vtkm::Id3(this->BlockSize - 1)) / vtkm::Id3(this->BlockSize);

  vtkm::cont::CellSetStructured<3> newCellSet;
  newCellSet.SetPointDimensions(numBlocks);

  if (!input.GetCoordinateSystem().GetData().CanConvert<vtkm::cont::ArrayHandleUniformPointCoordinates>())
  {
    // Technically, the filter would still "work", but the new coordinates below is only
    // valid with uniform coordiantes.
    throw vtkm::cont::ErrorBadType("Extraction only works with uniform point coordinates.");
  }
  vtkm::Bounds bounds = input.GetCoordinateSystem().GetBounds();
  vtkm::Vec3f origin{ bounds.MinCorner() };
  vtkm::Vec3f spacing{ (bounds.MaxCorner() - bounds.MinCorner()) / (numBlocks - 1) };
  vtkm::cont::ArrayHandleUniformPointCoordinates newCoordinates{ numBlocks, origin, spacing };

  // Create key that groups subsampling
  vtkm::cont::ArrayHandle<vtkm::Id> keyArray;
  this->Invoke(CreateNewKeyWorklet{numPoints, numBlocks, this->BlockSize},
               vtkm::cont::ArrayHandleIndex{ numPoints[0] * numPoints[1] * numPoints[2] },
               keyArray);

  vtkm::worklet::Keys<vtkm::Id> keys{ keyArray };
  auto mapper = [&](vtkm::cont::DataSet& data, const vtkm::cont::Field& field) {
    ComputeMeanStdevForField(this, data, field, keys);
  };
  return this->CreateResultCoordinateSystem(input,
                                            newCellSet,
                                            input.GetCoordinateSystem().GetName(),
                                            newCoordinates,
                                            mapper);
}

}
}
} // namespace vtkm::filter::uncertainty
