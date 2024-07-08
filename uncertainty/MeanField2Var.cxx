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
#include "./MeanField2Var.h"
#include "./worklet/MeanField2Var.h"

namespace vtkm
{
  namespace filter
  {
    namespace uncertainty
    {
      VTKM_CONT vtkm::cont::DataSet FiberMean::DoExecute(const vtkm::cont::DataSet &input)
      {

        // vtkm::cont::Timer timer;
        // std::cout << "detailed timer device: " << timer.GetDevice().GetName() << std::endl;
        // timer.Start();
        //  Input Field
        vtkm::cont::Field EnsembleMinOne = this->GetFieldFromDataSet(0, input);
        vtkm::cont::Field EnsembleMaxOne = this->GetFieldFromDataSet(1, input);
        vtkm::cont::Field EnsembleMinTwo = this->GetFieldFromDataSet(2, input);
        vtkm::cont::Field EnsembleMaxTwo = this->GetFieldFromDataSet(3, input);

        // Output Field
        vtkm::cont::UnknownArrayHandle OutputProbability;

        // timer.Stop();
        // std::cout << "filter 1 " << timer.GetElapsedTime() << std::endl;
        //  For Invoker
        auto resolveType = [&](auto ConcreteEnsembleMinOne)
        {
          // timer.Start();
          //  Obtaining Type
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
          vtkm::cont::ArrayHandle<ValueType> Probability;
          // timer.Stop();
          // std::cout << "filter 2 " << timer.GetElapsedTime() << std::endl;
          // timer.Start();

          // Invoker

          if (this->Approach == "MonteCarlo")
          {
            std::cout << "Adopt monte carlo with numsamples " << this->NumSamples << std::endl;
            this->Invoke(vtkm::worklet::detail::MultiVariateMonteCarlo{this->minAxis, this->maxAxis, this->NumSamples},
                         ConcreteEnsembleMinOne,
                         ConcreteEnsembleMaxOne,
                         ConcreteEnsembleMinTwo,
                         ConcreteEnsembleMaxTwo,
                         Probability);
          }
          else if (this->Approach == "ClosedForm")
          {
            std::cout << "Adopt ClosedForm" << std::endl;
            this->Invoke(vtkm::worklet::detail::MultiVariateClosedForm{this->minAxis, this->maxAxis},
                         ConcreteEnsembleMinOne,
                         ConcreteEnsembleMaxOne,
                         ConcreteEnsembleMinTwo,
                         ConcreteEnsembleMaxTwo,
                         Probability);
          }
          else if (this->Approach == "Mean")
          {
            std::cout << "Adopt Mean" << std::endl;
            this->Invoke(vtkm::worklet::detail::MultiVariateMean{this->minAxis, this->maxAxis},
                         ConcreteEnsembleMinOne,
                         ConcreteEnsembleMaxOne,
                         ConcreteEnsembleMinTwo,
                         ConcreteEnsembleMaxTwo,
                         Probability);
          }
          else
          {
            throw std::runtime_error("unsupported approach:" + this->Approach);
          }

          // From Temporary Output Variable to Output Variable
          OutputProbability = Probability;
          // timer.Stop();
          // std::cout << "filter 3 " << timer.GetElapsedTime() << std::endl;
        };
        this->CastAndCallScalarField(EnsembleMinOne, resolveType);

        // Creating Result
        // timer.Start();
        vtkm::cont::DataSet result = this->CreateResult(input);
        result.AddPointField("OutputProbability", OutputProbability);
        // timer.Stop();
        // std::cout << "filter 4 " << timer.GetElapsedTime() << std::endl;

        return result;
      }
    }
  }
}
