//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include "ContourUncertainEnsemble2D.h"

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>

// #include "ucvworklet/MVGaussianWithEnsemble3D.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2DTryLialgEntropy.hpp"

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

//TODO, more flexible way to support different legnth of vector
using SupportedTypesVec = vtkm::List<vtkm::Vec<double, 20>,vtkm::Vec<double, 15>>;

ContourUncertainEnsemble2D::ContourUncertainEnsemble2D()
{
  this->SetCrossProbabilityName("cross_probability");
}

vtkm::cont::DataSet ContourUncertainEnsemble2D::DoExecute(const vtkm::cont::DataSet& input)
{
  
  vtkm::cont::Field ensembleField = this->GetFieldFromDataSet(0, input);
  //vtkm::cont::Field meanField = this->GetFieldFromDataSet(1, input);
  
  vtkm::cont::UnknownArrayHandle crossProbability;
  vtkm::cont::UnknownArrayHandle numNonZeroProbability;
  vtkm::cont::UnknownArrayHandle entropy;

  if (!input.GetCellSet().IsType<vtkm::cont::CellSetStructured<2>>())
  {
    input.GetCellSet().PrintSummary(std::cout);
    throw vtkm::cont::ErrorBadType("Uncertain contour only works for CellSetStructured data now.");
  }

  auto resolveType = [&](auto concrete) {

  vtkm::cont::ArrayHandle<vtkm::FloatDefault> concreteCrossProb;
  vtkm::cont::ArrayHandle<vtkm::Id> concreteNumNonZeroProb;
  vtkm::cont::ArrayHandle<vtkm::FloatDefault> concreteEntropy;

  this->Invoke(MVGaussianWithEnsemble2DTryLialgEntropy{ this->IsoValue, 1000 },
                 input.GetCellSet(),
                 concrete,
                 concreteCrossProb,
                 concreteNumNonZeroProb,
                 concreteEntropy);

    crossProbability = concreteCrossProb;
    numNonZeroProbability = concreteNumNonZeroProb;
    entropy = concreteEntropy;
  };
  //this->CastAndCallScalarField(ensembleField, resolveType);
  ensembleField.GetData().CastAndCallForTypes<SupportedTypesVec, VTKM_DEFAULT_STORAGE_LIST>(resolveType);
  vtkm::cont::DataSet result = this->CreateResult(input);
  result.AddCellField(this->GetCrossProbabilityName(), crossProbability);
  result.AddCellField(this->GetNumberNonzeroProbabilityName(), numNonZeroProbability);
  result.AddCellField(this->GetEntropyName(), entropy);
  return result;
}

}
}
} // vtkm::filter::uncertainty
