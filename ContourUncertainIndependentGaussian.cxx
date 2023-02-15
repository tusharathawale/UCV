//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "ContourUncertainIndependentGaussian.h"

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>

#include "ucvworklet/EntropyIndependentGaussian.hpp"

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

ContourUncertainIndependentGaussian::ContourUncertainIndependentGaussian()
{
  this->SetCrossProbabilityName("cross_probability");
}

vtkm::cont::DataSet ContourUncertainIndependentGaussian::DoExecute(const vtkm::cont::DataSet& input)
{
  vtkm::cont::Field meanField = this->GetFieldFromDataSet(0, input);
  vtkm::cont::Field stdevField = this->GetFieldFromDataSet(1, input);

  vtkm::cont::UnknownArrayHandle crossProbability;
  vtkm::cont::UnknownArrayHandle numNonZeroProbability;
  vtkm::cont::UnknownArrayHandle entropy;

  if (!input.GetCellSet().IsType<vtkm::cont::CellSetStructured<3>>())
  {
    throw vtkm::cont::ErrorBadType("Uncertain contour only works for CellSetStructured<3>.");
  }
  vtkm::cont::CellSetStructured<3> cellSet;
  input.GetCellSet().AsCellSet(cellSet);

  auto resolveType = [&](auto concreteMeanField) {
    using ArrayType = std::decay_t<decltype(concreteMeanField)>;
    using ValueType = typename ArrayType::ValueType;
    ArrayType concreteStdevField;
    vtkm::cont::ArrayCopyShallowIfPossible(stdevField.GetData(), concreteStdevField);

    vtkm::cont::ArrayHandle<ValueType> concreteCrossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> concreteNumNonZeroProb;
    vtkm::cont::ArrayHandle<ValueType> concreteEntropy;

    this->Invoke(EntropyIndependentGaussian{ this->IsoValue },
                 cellSet,
                 concreteMeanField,
                 concreteStdevField,
                 concreteCrossProb,
                 concreteNumNonZeroProb,
                 concreteEntropy);

    crossProbability = concreteCrossProb;
    numNonZeroProbability = concreteNumNonZeroProb;
    entropy = concreteEntropy;
  };
  this->CastAndCallScalarField(meanField, resolveType);

  vtkm::cont::DataSet result = this->CreateResult(input);
  result.AddCellField(this->GetCrossProbabilityName(), crossProbability);
  result.AddCellField(this->GetNumberNonzeroProbabilityName(), numNonZeroProbability);
  result.AddCellField(this->GetEntropyName(), entropy);
  return result;
}

}
}
} // vtkm::filter::uncertainty
