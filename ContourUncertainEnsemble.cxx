//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "ContourUncertainEnsemble.h"

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>

// #include "ucvworklet/MVGaussianWithEnsemble3D.hpp"
#include "ucvworklet/MVGaussianWithEnsemble3DTryLialg.hpp"

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

ContourUncertainEnsemble::ContourUncertainEnsemble()
{
  this->SetCrossProbabilityName("cross_probability");
}

vtkm::cont::DataSet ContourUncertainEnsemble::DoExecute(const vtkm::cont::DataSet& input)
{
  vtkm::cont::Field ensembleField = this->GetFieldFromDataSet(0, input);
  vtkm::cont::Field meanField = this->GetFieldFromDataSet(1, input);

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
    using ValueType = typename std::decay_t<decltype(concreteMeanField)>::ValueType;
    vtkm::cont::ArrayHandle<vtkm::Vec<ValueType, 4 * 4 * 4>> concreteEnsembleField;
    vtkm::cont::ArrayCopyShallowIfPossible(ensembleField.GetData(), concreteEnsembleField);

    vtkm::cont::ArrayHandle<ValueType> concreteCrossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> concreteNumNonZeroProb;
    vtkm::cont::ArrayHandle<ValueType> concreteEntropy;

    this->Invoke(MVGaussianWithEnsemble3DTryLialg{ this->IsoValue, 1000 },
                 cellSet,
                 concreteEnsembleField,
                 concreteMeanField,
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
