//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "SubsampleUncertaintyEnsemble.h"

#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayHandleUniformPointCoordinates.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ErrorBadType.h>
#include <vtkm/cont/ErrorBadValue.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/Keys.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/ExtractingMeanRaw.hpp"

constexpr vtkm::IdComponent FORCE_BLOCK_SIZE = 4;
constexpr vtkm::IdComponent FORCE_ENSEMBLE_SIZE =
  FORCE_BLOCK_SIZE * FORCE_BLOCK_SIZE * FORCE_BLOCK_SIZE;

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

vtkm::cont::DataSet SubsampleUncertaintyEnsemble::DoExecute(const vtkm::cont::DataSet& input)
{
  if (!input.GetCellSet().IsType<vtkm::cont::CellSetStructured<3>>())
  {
    throw vtkm::cont::ErrorBadType("Extraction only works for CellSetStructured<3>.");
  }
  vtkm::cont::CellSetStructured<3> cellSet;
  input.GetCellSet().AsCellSet(cellSet);

  vtkm::Id3 numPoints = cellSet.GetPointDimensions();

  vtkm::Id3 numBlocks = (numPoints + vtkm::Id3(this->BlockSize - 1)) / vtkm::Id3(this->BlockSize);

  // Currently enforce a blocksize of 4 and a cell set with dimensions evenly divisible by the
  // blocksize. We could support general Ensebles by putting them in an ArrayHandleGroupVecVariable,
  // but that is not supported by ContourUncertainEnsemble, which is currently the only point
  // of this code.
  if (this->BlockSize != FORCE_BLOCK_SIZE)
  {
    throw vtkm::cont::ErrorBadValue("Multivariant Gaussian currently only supports block size of 4.");
  }
  if (((numPoints[0] % 4) != 0) || ((numPoints[1] % 4) != 0) || ((numPoints[2] % 4) != 0))
  {
    throw vtkm::cont::ErrorBadValue("Multivariant Gaussian currently only supports grids divisible by 4 in each dimension.");
  }

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
    this->MapField(data, field, keys);
  };
  return this->CreateResultCoordinateSystem(input,
                                            newCellSet,
                                            input.GetCoordinateSystem().GetName(),
                                            newCoordinates,
                                            mapper);
}

VTKM_CONT void SubsampleUncertaintyEnsemble::MapField(
    vtkm::cont::DataSet& data,
    const vtkm::cont::Field& field,
    const vtkm::worklet::Keys<vtkm::Id>& keys) const
{
  if (field.IsPointField())
  {
    vtkm::cont::UnknownArrayHandle meanArray;
    vtkm::cont::UnknownArrayHandle ensembleArray;
    // This code only supports scalar fields and converts fields to floating point. It would
    // be possible to preserve the original data type, but ContourUncertainEnsemble only works
    // with float arrays. Plus, the mean might be less accurate with integer types. It would
    // also be possible to support arbitrary vectors (especially with some new features in
    // VTK-m 2.1), but, again, this is not supported by ContourUncertainEnsemble, so we don't.
    auto resolveType = [&](const auto& concrete) {
      using ValueType = typename std::decay_t<decltype(concrete)>::ValueType;
      vtkm::cont::ArrayHandle<ValueType> meanConcrete;
      vtkm::cont::ArrayHandle<vtkm::Vec<ValueType, FORCE_ENSEMBLE_SIZE>> ensembleConcrete;
      this->Invoke(ExtractingMeanRaw{}, keys, concrete, meanConcrete, ensembleConcrete);
      meanArray = meanConcrete;
      ensembleArray = ensembleConcrete;
    };
    // This would be easier with CastAndCallScalarField, but this is only available in
    // FilterField. Then again, the use in MapField is not great as it is usually better
    // to preserve the original data type.
    field.GetData().CastAndCallForTypesWithFloatFallback<vtkm::TypeListFieldScalar, VTKM_DEFAULT_STORAGE_LIST>(
          resolveType);
    data.AddPointField(field.GetName(), meanArray);
    data.AddPointField(field.GetName() + this->GetEnsembleSuffix(), ensembleArray);
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

}
}
} // namespace vtkm::filter::uncertainty
