//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "ContourUncertainUniform.h"

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>

#include "ucvworklet/EntropyUniform.hpp"

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

ContourUncertainUniform::ContourUncertainUniform()
{
  this->SetCrossProbabilityName("cross_probability");
}

vtkm::cont::DataSet ContourUncertainUniform::DoExecute(const vtkm::cont::DataSet& input)
{
  vtkm::cont::Field minField = this->GetFieldFromDataSet(0, input);
  vtkm::cont::Field maxField = this->GetFieldFromDataSet(1, input);

  vtkm::cont::UnknownArrayHandle crossProbability;
  vtkm::cont::UnknownArrayHandle numNonZeroProbability;
  vtkm::cont::UnknownArrayHandle entropy;

  if (!input.GetCellSet().IsType<vtkm::cont::CellSetStructured<3>>())
  {
    throw vtkm::cont::ErrorBadType("Uncertain contour only works for CellSetStructured<3>.");
  }
  vtkm::cont::CellSetStructured<3> cellSet;
  input.GetCellSet().AsCellSet(cellSet);

  auto resolveType = [&](auto concreteMinField) {
    using ArrayType = std::decay_t<decltype(concreteMinField)>;
    using ValueType = typename ArrayType::ValueType;
    ArrayType concreteMaxField;
    vtkm::cont::ArrayCopyShallowIfPossible(maxField.GetData(), concreteMaxField);

    vtkm::cont::ArrayHandle<ValueType> concreteCrossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> concreteNumNonZeroProb;
    vtkm::cont::ArrayHandle<ValueType> concreteEntropy;

    this->Invoke(EntropyUniform{ this->IsoValue },
                 cellSet,
                 concreteMinField,
                 concreteMaxField,
                 concreteCrossProb,
                 concreteNumNonZeroProb,
                 concreteEntropy);

    crossProbability = concreteCrossProb;
    numNonZeroProbability = concreteNumNonZeroProb;
    entropy = concreteEntropy;
  };
  this->CastAndCallScalarField(minField, resolveType);

  vtkm::cont::DataSet result = this->CreateResult(input);
  result.AddCellField(this->GetCrossProbabilityName(), crossProbability);
  result.AddCellField(this->GetNumberNonzeroProbabilityName(), numNonZeroProbability);
  result.AddCellField(this->GetEntropyName(), entropy);
  return result;
}

}
}
} // vtkm::filter::uncertainty
