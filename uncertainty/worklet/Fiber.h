//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
// New Fiber.h

#ifndef vtk_m_worklet_uncertainty_Fiber_inner_h
#define vtk_m_worklet_uncertainty_Fiber_inner_h
#include <iostream>
#include <utility>
#include <vector>
#include <vtkm/worklet/WorkletMapField.h>

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
#include <thrust/device_vector.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/uniform_real_distribution.h>
#else
#include <random>
#endif

namespace vtkm
{
  namespace worklet
  {
    namespace detail
    {
      class FiberMonteCarlo : public vtkm::worklet::WorkletMapField
      {
      public:
        // Worklet Input
        // Fiber(const std::vector<std::pair<double, double>>& minAxis,
        //      const std::vector<std::pair<double, double>>& maxAxis)
        //  : InputMinAxis(minAxis), InputMaxAxis(maxAxis){};
        FiberMonteCarlo(const vtkm::Pair<vtkm::Float64, vtkm::Float64> &minAxis,
                        const vtkm::Pair<vtkm::Float64, vtkm::Float64> &maxAxis,
                        const vtkm::Id numSamples)
            : InputMinAxis(minAxis), InputMaxAxis(maxAxis), NumSamples(numSamples){};

        // Input and Output Parameters
        using ControlSignature = void(FieldIn, FieldIn, FieldIn, FieldIn, FieldOut);

        using ExecutionSignature = void(_1, _2, _3, _4, _5);
        // using ExecutionSignature = void(_2, _3, _4, _5);
        using InputDomain = _1;

        // Template
        template <typename MinOne,
                  typename MaxOne,
                  typename MinTwo,
                  typename MaxTwo,
                  typename OutCellFieldType>
        // Operator
        VTKM_EXEC void operator()(const MinOne &EnsembleMinOne,
                                  const MaxOne &EnsembleMaxOne,
                                  const MinTwo &EnsembleMinTwo,
                                  const MaxTwo &EnsembleMaxTwo,
                                  OutCellFieldType &probability) const
        {
          // User defined rectangle(trait)
          vtkm::FloatDefault X1 = 0.0;
          X1 = static_cast<vtkm::FloatDefault>(InputMinAxis.first);
          vtkm::FloatDefault Y1 = 0.0;
          Y1 = static_cast<vtkm::FloatDefault>(InputMinAxis.second);
          vtkm::FloatDefault X2 = 0.0;
          X2 = static_cast<vtkm::FloatDefault>(InputMaxAxis.first);
          vtkm::FloatDefault Y2 = 0.0;
          Y2 = static_cast<vtkm::FloatDefault>(InputMaxAxis.second);
          vtkm::FloatDefault TraitArea = (X2 - X1) * (Y2 - Y1);

          vtkm::FloatDefault X5 = 0.0;
          vtkm::FloatDefault X6 = 0.0;
          vtkm::FloatDefault Y5 = 0.0;
          vtkm::FloatDefault Y6 = 0.0;

          vtkm::FloatDefault X3 = 0.0;
          vtkm::FloatDefault Y3 = 0.0;
          vtkm::FloatDefault X4 = 0.0;
          vtkm::FloatDefault Y4 = 0.0;

          vtkm::FloatDefault IntersectionArea = 0.0;
          vtkm::FloatDefault IntersectionProbablity = 0.0;
          vtkm::FloatDefault IntersectionHeight = 0.0;
          vtkm::FloatDefault IntersectionWidth = 0.0;

          vtkm::FloatDefault N1 = 0.0;
          vtkm::FloatDefault N2 = 0.0;
          vtkm::IdComponent NonZeroCases = 0;
          vtkm::FloatDefault MCProbability = 0.0;

          // data rectangle
          X3 = static_cast<vtkm::FloatDefault>(EnsembleMinOne);
          X4 = static_cast<vtkm::FloatDefault>(EnsembleMaxOne);
          Y3 = static_cast<vtkm::FloatDefault>(EnsembleMinTwo);
          Y4 = static_cast<vtkm::FloatDefault>(EnsembleMaxTwo);

          // Monte Carlo
          // Trait Coordinates (X1,Y1) & (X2,Y2)
#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
          thrust::minstd_rand rng;
          thrust::uniform_real_distribution<vtkm::FloatDefault> distX(X1, X2);
          thrust::uniform_real_distribution<vtkm::FloatDefault> distY(Y1, Y2);

          for (vtkm::IdComponent i = 0; i < NumSample; i++)
          {
            N1 = distX(rng);
            N2 = distY(rng);
            if ((N1 > X3) && (N1 < X4) && (N2 > Y3) && (N2 < Y4))
            {
              NonZeroCases++;
            }
          }

#else
          std::random_device rd;
          std::mt19937 gen(rd());
          // Generate samples from data rectangle
          std::uniform_real_distribution<vtkm::FloatDefault> GenerateN1(X3, X4);
          std::uniform_real_distribution<vtkm::FloatDefault> GenerateN2(Y3, Y4);

          for (vtkm::IdComponent i = 0; i < this->NumSamples; i++)
          {
            N1 = GenerateN1(gen);
            N2 = GenerateN2(gen);
            // increase the case number when the data is located in user rectangle
            if ((N1 > X1) and (N1 < X2) and (N2 > Y1) and (N2 < Y2))
            {
              NonZeroCases++;
            }
          }

#endif
          //printf("NonZeroCases %d this->NumSamples %d\n", NonZeroCases, this->NumSamples);
          // printf("X1 %f Y1 %f X2 %f Y2 %f X3 %f X4 %f Y3 %f Y4 %f NonZeroCases %d\n",X1,Y1,X2,Y2,X3,X4,Y3,Y4,NonZeroCases);
          MCProbability = 1.0 * NonZeroCases / (1.0 * this->NumSamples);
          probability = MCProbability;

          return;
        }

      private:
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputMinAxis;
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputMaxAxis;
        vtkm::Id NumSamples = 1;
      };

      class FiberClosedForm : public vtkm::worklet::WorkletMapField
      {
      public:
        // Worklet Input
        // Fiber(const std::vector<std::pair<double, double>>& minAxis,
        //      const std::vector<std::pair<double, double>>& maxAxis)
        //  : InputMinAxis(minAxis), InputMaxAxis(maxAxis){};
        FiberClosedForm(const vtkm::Pair<vtkm::Float64, vtkm::Float64> &minAxis,
                        const vtkm::Pair<vtkm::Float64, vtkm::Float64> &maxAxis)
            : InputMinAxis(minAxis), InputMaxAxis(maxAxis){};

        // Input and Output Parameters
        using ControlSignature = void(FieldIn, FieldIn, FieldIn, FieldIn, FieldOut);

        using ExecutionSignature = void(_1, _2, _3, _4, _5);
        // using ExecutionSignature = void(_2, _3, _4, _5);
        using InputDomain = _1;

        // Template
        template <typename MinOne,
                  typename MaxOne,
                  typename MinTwo,
                  typename MaxTwo,
                  typename OutCellFieldType>
        // Operator
        VTKM_EXEC void operator()(const MinOne &EnsembleMinOne,
                                  const MaxOne &EnsembleMaxOne,
                                  const MinTwo &EnsembleMinTwo,
                                  const MaxTwo &EnsembleMaxTwo,
                                  OutCellFieldType &probability) const
        {
          // User defined rectangle(trait)
          vtkm::FloatDefault X1 = 0.0;
          X1 = static_cast<vtkm::FloatDefault>(InputMinAxis.first);
          vtkm::FloatDefault Y1 = 0.0;
          Y1 = static_cast<vtkm::FloatDefault>(InputMinAxis.second);
          vtkm::FloatDefault X2 = 0.0;
          X2 = static_cast<vtkm::FloatDefault>(InputMaxAxis.first);
          vtkm::FloatDefault Y2 = 0.0;
          Y2 = static_cast<vtkm::FloatDefault>(InputMaxAxis.second);
          vtkm::FloatDefault TraitArea = (X2 - X1) * (Y2 - Y1);

          vtkm::FloatDefault X5 = 0.0;
          vtkm::FloatDefault X6 = 0.0;
          vtkm::FloatDefault Y5 = 0.0;
          vtkm::FloatDefault Y6 = 0.0;

          vtkm::FloatDefault X3 = 0.0;
          vtkm::FloatDefault Y3 = 0.0;
          vtkm::FloatDefault X4 = 0.0;
          vtkm::FloatDefault Y4 = 0.0;

          vtkm::FloatDefault IntersectionArea = 0.0;
          vtkm::FloatDefault IntersectionProbablity = 0.0;
          vtkm::FloatDefault IntersectionHeight = 0.0;
          vtkm::FloatDefault IntersectionWidth = 0.0;

          vtkm::FloatDefault N1 = 0.0;
          vtkm::FloatDefault N2 = 0.0;
          vtkm::IdComponent NonZeroCases = 0;
          vtkm::IdComponent NumSample = 500;
          vtkm::FloatDefault MCProbability = 0.0;

          // data rectangle
          X3 = static_cast<vtkm::FloatDefault>(EnsembleMinOne);
          X4 = static_cast<vtkm::FloatDefault>(EnsembleMaxOne);
          Y3 = static_cast<vtkm::FloatDefault>(EnsembleMinTwo);
          Y4 = static_cast<vtkm::FloatDefault>(EnsembleMaxTwo);

          // caculating intersection of two rectangle (overlapping region)
          X5 = std::max(X1, X3);
          Y5 = std::max(Y1, Y3);
          X6 = std::min(X2, X4);
          Y6 = std::min(Y2, Y4);

          IntersectionHeight = Y6 - Y5;
          IntersectionWidth = X6 - X5;

          vtkm::FloatDefault DataArea = (X4 - X3) * (Y4 - Y3);

          if ((IntersectionHeight > 0) and (IntersectionWidth > 0) and (X5 < X6) and (Y5 < Y6))
          {
            IntersectionArea = IntersectionHeight * IntersectionWidth;
            // the portion of trait
            // IntersectionProbablity = IntersectionArea / TraitArea;
            IntersectionProbablity = IntersectionArea / DataArea;
          }
          probability = IntersectionProbablity;
          return;
        }

      private:
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputMinAxis;
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputMaxAxis;
      };
    }

  }
}
#endif
