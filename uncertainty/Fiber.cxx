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
#include "./Fiber.h"
#include "./worklet/Fiber.h"

namespace vtkm
{
  namespace filter
  {
    namespace uncertainty
    {
      VTKM_CONT vtkm::cont::DataSet Fiber::DoExecute(const vtkm::cont::DataSet &input)
      {

        vtkm::cont::Timer timer;
        std::cout << "detailed timer device: " << timer.GetDevice().GetName() << std::endl;
        timer.Start();
        // Input Field
        vtkm::cont::Field EnsembleMinOne = this->GetFieldFromDataSet(0, input);
        vtkm::cont::Field EnsembleMaxOne = this->GetFieldFromDataSet(1, input);
        vtkm::cont::Field EnsembleMinTwo = this->GetFieldFromDataSet(2, input);
        vtkm::cont::Field EnsembleMaxTwo = this->GetFieldFromDataSet(3, input);

        // Output Field
        vtkm::cont::UnknownArrayHandle OutputMonteCarloProbability;
        vtkm::cont::UnknownArrayHandle OutputInteriorProbability;

        timer.Stop();
        std::cout << "filter 1 " << timer.GetElapsedTime() << std::endl;
        // For Invoker
        auto resolveType = [&](auto ConcreteEnsembleMinOne)
        {
          timer.Start();
          // Obtaining Type
          using ArrayType = std::decay_t<decltype(ConcreteEnsembleMinOne)>;
          using ValueType = typename ArrayType::ValueType;

          // Temporary Input Variable to add input values
          ArrayType ConcreteEnsembleMaxOne;
          ArrayType ConcreteEnsembleMinTwo;
          ArrayType ConcreteEnsembleMaxTwo;

          vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMaxOne.GetData(), ConcreteEnsembleMaxOne);
          vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMinTwo.GetData(), ConcreteEnsembleMinTwo);
          vtkm::cont::ArrayCopyShallowIfPossible(EnsembleMaxTwo.GetData(), ConcreteEnsembleMaxTwo);

          // Temporary Output Variable
          vtkm::cont::ArrayHandle<ValueType> ConcreteMonteCarloProbability;
          vtkm::cont::ArrayHandle<ValueType> ConcreteInteriorProbability;
          timer.Stop();
          std::cout << "filter 2 " << timer.GetElapsedTime() << std::endl;
          timer.Start();

          // Invoker
          // this->IsoValue
          this->Invoke(vtkm::worklet::detail::Fiber{this->minAxis, this->maxAxis},
                       ConcreteEnsembleMinOne,
                       ConcreteEnsembleMaxOne,
                       ConcreteEnsembleMinTwo,
                       ConcreteEnsembleMaxTwo,
                       ConcreteMonteCarloProbability,
                       ConcreteInteriorProbability);

          // From Temporary Output Variable to Output Variable
          OutputMonteCarloProbability = ConcreteMonteCarloProbability;
          OutputInteriorProbability = ConcreteInteriorProbability;
          timer.Stop();
          std::cout << "filter 3 " << timer.GetElapsedTime() << std::endl;
        };
        this->CastAndCallScalarField(EnsembleMinOne, resolveType);
        
        // Creating Result
        timer.Start();
        vtkm::cont::DataSet result = this->CreateResult(input);
        result.AddPointField("OutputMonteCarloProbability", OutputMonteCarloProbability);
        result.AddPointField("OutputInteriorProbability", OutputInteriorProbability);
        timer.Stop();
        std::cout << "filter 4 " << timer.GetElapsedTime() << std::endl;

        return result;
      }
    }
  }
}
