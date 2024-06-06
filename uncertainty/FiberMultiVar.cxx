//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/filter/uncertainty/FiberMultiVar.h>
#include <vtkm/filter/uncertainty/worklet/FiberMultiVar.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{
VTKM_CONT vtkm::cont::DataSet FiberMultiVar::DoExecute(const vtkm::cont::DataSet& input)
{
  // vtkm::cont::Timer timer;
  // timer.Start();

  // Input Field
  vtkm::cont::Field EnsembleMinX = this->GetFieldFromDataSet(0, input); // For resolve type
  vtkm::cont::Field EnsembleMaxX = this->GetFieldFromDataSet(1, input);
  vtkm::cont::Field EnsembleMinY = this->GetFieldFromDataSet(2, input);
  vtkm::cont::Field EnsembleMaxY = this->GetFieldFromDataSet(3, input);
  vtkm::cont::Field EnsembleMinZ = this->GetFieldFromDataSet(4, input); 
  vtkm::cont::Field EnsembleMaxZ = this->GetFieldFromDataSet(5, input);


  // Output Field
  vtkm::cont::UnknownArrayHandle OutputMonteCarloProbability;
  vtkm::cont::UnknownArrayHandle OutputInteriorProbability;
  
  // CellSet
  vtkm::cont::CellSetStructured<3> cellSet;
  input.GetCellSet().AsCellSet(cellSet);


  // For Invoker
  auto resolveType = [&](auto ConcreteEnsembleMinX) {
    
    // Obtaining Type
    using ArrayType = std::decay_t<decltype(ConcreteEnsembleMinX)>;
    using ValueType = typename ArrayType::ValueType;

    // Temporary Input Variable to add input values
    ArrayType ConcreteEnsembleMaxX;
    ArrayType ConcreteEnsembleMinY;
    ArrayType ConcreteEnsembleMaxY;
    ArrayType ConcreteEnsembleMinZ;
    ArrayType ConcreteEnsembleMaxZ;

    vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMaxX.GetData(), ConcreteEnsembleMaxX);
    vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMinY.GetData(), ConcreteEnsembleMinY);
    vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMaxY.GetData(), ConcreteEnsembleMaxY);
    vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMinZ.GetData(), ConcreteEnsembleMinZ);
    vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMaxZ.GetData(), ConcreteEnsembleMaxZ);


    // Temporary Output Variable
    vtkm::cont::ArrayHandle<ValueType> ConcreteMonteCarloProbability;
    vtkm::cont::ArrayHandle<ValueType> ConcreteInteriorProbability;

    // Invoker
    // this->IsoValue
    this->Invoke(vtkm::worklet::detail::FiberMultiVar{ this->bottomLeft, this->topRight },
                 cellSet,
                 ConcreteEnsembleMinX,
                 ConcreteEnsembleMaxX,
                 ConcreteEnsembleMinY,
                 ConcreteEnsembleMaxY,
                 ConcreteEnsembleMinZ,
                 ConcreteEnsembleMaxZ,
                 ConcreteMonteCarloProbability,
                 ConcreteInteriorProbability);

    // From Temporary Output Variable to Output Variable
    OutputMonteCarloProbability = ConcreteMonteCarloProbability;
    OutputInteriorProbability = ConcreteInteriorProbability;
  };
  this->CastAndCallScalarField(EnsembleMinX, resolveType);
  
  // Creating Result
  vtkm::cont::DataSet result = this->CreateResult(input);
  result.AddPointField("OutputMonteCarloProbability", OutputMonteCarloProbability);
  result.AddPointField("OutputInteriorProbability", OutputInteriorProbability);
  
  //   timer.Stop();
  //   vtkm::Float64 elapsedTime = timer.GetElapsedTime();
  //   std::cout << elapsedTime << std::endl;
  
  return result;
}
}
}
}
