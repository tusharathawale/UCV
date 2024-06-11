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
        //  : InputBottomLeft(minAxis), InputTopRight(maxAxis){};
        FiberMonteCarlo(const vtkm::Pair<vtkm::Float64, vtkm::Float64> &minAxis,
                        const vtkm::Pair<vtkm::Float64, vtkm::Float64> &maxAxis,
                        const vtkm::Id numSamples)
            : InputBottomLeft(minAxis), InputTopRight(maxAxis), NumSamples(numSamples){};

        // Input and Output Parameters
        using ControlSignature = void(FieldIn, FieldIn, FieldIn, FieldIn, FieldOut);

        using ExecutionSignature = void(_1, _2, _3, _4, _5);
        // using ExecutionSignature = void(_2, _3, _4, _5);
        using InputDomain = _1;

        // Template
        template <typename MinX,
                  typename MaxX,
                  typename MinY,
                  typename MaxY,
                  typename OutCellFieldType>
        // Operator
        VTKM_EXEC void operator()(const MinX &EnsembleMinX,
                                  const MaxX &EnsembleMaxX,
                                  const MinY &EnsembleMinY,
                                  const MaxY &EnsembleMaxY,
                                  OutCellFieldType &probability) const
        {
          // User defined rectangle(trait)
          vtkm::FloatDefault minX_user = 0.0;
          minX_user = static_cast<vtkm::FloatDefault>(InputBottomLeft.first);
          vtkm::FloatDefault minY_user = 0.0;
          minY_user = static_cast<vtkm::FloatDefault>(InputBottomLeft.second);
          vtkm::FloatDefault maxX_user = 0.0;
          maxX_user = static_cast<vtkm::FloatDefault>(InputTopRight.first);
          vtkm::FloatDefault maxY_user = 0.0;
          maxY_user = static_cast<vtkm::FloatDefault>(InputTopRight.second);
          // vtkm::FloatDefault TraitArea = (maxX_user - minX_user) * (maxY_user - minY_user);

          vtkm::FloatDefault minX_dataset = 0.0;
          vtkm::FloatDefault minY_dataset = 0.0;
          vtkm::FloatDefault maxX_dataset = 0.0;
          vtkm::FloatDefault maxY_dataset = 0.0;

          vtkm::FloatDefault N1 = 0.0;
          vtkm::FloatDefault N2 = 0.0;
          vtkm::IdComponent NonZeroCases = 0;
          vtkm::FloatDefault MCProbability = 0.0;

          // data rectangle
          minX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinX);
          maxX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxX);
          minY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinY);
          maxY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxY);

          // if data rectangle is zero, there is no uncertainty, return zero
          if (abs(minX_dataset - maxX_dataset) < 0.000001 && abs(minY_dataset - maxY_dataset) < 0.000001)
          {
            probability = 0.0;
            return;
          }

          // Monte Carlo
          // Trait Coordinates (minX_user,minY_user) & (maxX_user,maxY_user)
#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
          thrust::minstd_rand rng;
          thrust::uniform_real_distribution<vtkm::FloatDefault> distX(minX_dataset, maxX_dataset);
          thrust::uniform_real_distribution<vtkm::FloatDefault> distY(minY_dataset, maxY_dataset);

          for (vtkm::IdComponent i = 0; i < this->NumSamples; i++)
          {
            N1 = distX(rng);
            N2 = distY(rng);
            if ((N1 > minX_user) && (N1 < maxX_user) && (N2 > minY_user) && (N2 < maxY_user))
            {
              NonZeroCases++;
            }
          }

#else
          std::random_device rd;
          std::mt19937 gen(rd());
          // Generate samples from data rectangle
          std::uniform_real_distribution<vtkm::FloatDefault> GenerateN1(minX_dataset, maxX_dataset);
          std::uniform_real_distribution<vtkm::FloatDefault> GenerateN2(minY_dataset, maxY_dataset);

          for (vtkm::IdComponent i = 0; i < this->NumSamples; i++)
          {
            N1 = GenerateN1(gen);
            N2 = GenerateN2(gen);
            // increase the case number when the data is located in user rectangle
            if ((N1 > minX_user) && (N1 < maxX_user) && (N2 > minY_user) && (N2 < maxY_user))
            {
              NonZeroCases++;
            }
          }

#endif
          // printf("NonZeroCases %d this->NumSamples %d\n", NonZeroCases, this->NumSamples);
          //  printf("minX_user %f minY_user %f maxX_user %f maxY_user %f minX_dataset %f maxX_dataset %f minY_dataset %f maxY_dataset %f NonZeroCases %d\n",minX_user,minY_user,maxX_user,maxY_user,minX_dataset,maxX_dataset,minY_dataset,maxY_dataset,NonZeroCases);
          MCProbability = 1.0 * NonZeroCases / (1.0 * this->NumSamples);
          probability = MCProbability;

          return;
        }

      private:
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputBottomLeft;
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputTopRight;
        vtkm::Id NumSamples = 1;
      };

      class FiberClosedForm : public vtkm::worklet::WorkletMapField
      {
      public:
        // Worklet Input
        // Fiber(const std::vector<std::pair<double, double>>& minAxis,
        //      const std::vector<std::pair<double, double>>& maxAxis)
        //  : InputBottomLeft(minAxis), InputTopRight(maxAxis){};
        FiberClosedForm(const vtkm::Pair<vtkm::Float64, vtkm::Float64> &minAxis,
                        const vtkm::Pair<vtkm::Float64, vtkm::Float64> &maxAxis)
            : InputBottomLeft(minAxis), InputTopRight(maxAxis){};

        // Input and Output Parameters
        using ControlSignature = void(FieldIn, FieldIn, FieldIn, FieldIn, FieldOut);

        using ExecutionSignature = void(_1, _2, _3, _4, _5);
        // using ExecutionSignature = void(_2, _3, _4, _5);
        using InputDomain = _1;

        // Template
        template <typename MinX,
                  typename MaxX,
                  typename MinY,
                  typename MaxY,
                  typename OutCellFieldType>
        // Operator
        VTKM_EXEC void operator()(const MinX &EnsembleMinX,
                                  const MaxX &EnsembleMaxX,
                                  const MinY &EnsembleMinY,
                                  const MaxY &EnsembleMaxY,
                                  OutCellFieldType &probability) const
        {
          // User defined rectangle(trait)
          vtkm::FloatDefault minX_user = 0.0;
          minX_user = static_cast<vtkm::FloatDefault>(InputBottomLeft.first);
          vtkm::FloatDefault minY_user = 0.0;
          minY_user = static_cast<vtkm::FloatDefault>(InputBottomLeft.second);
          vtkm::FloatDefault maxX_user = 0.0;
          maxX_user = static_cast<vtkm::FloatDefault>(InputTopRight.first);
          vtkm::FloatDefault maxY_user = 0.0;
          maxY_user = static_cast<vtkm::FloatDefault>(InputTopRight.second);
          // vtkm::FloatDefault TraitArea = (maxX_user - minX_user) * (maxY_user - minY_user);

          vtkm::FloatDefault minX_intersection = 0.0;
          vtkm::FloatDefault maxX_intersection = 0.0;
          vtkm::FloatDefault minY_intersection = 0.0;
          vtkm::FloatDefault maxY_intersection = 0.0;

          vtkm::FloatDefault minX_dataset = 0.0;
          vtkm::FloatDefault minY_dataset = 0.0;
          vtkm::FloatDefault maxX_dataset = 0.0;
          vtkm::FloatDefault maxY_dataset = 0.0;

          vtkm::FloatDefault IntersectionArea = 0.0;
          vtkm::FloatDefault IntersectionProbablity = 0.0;
          vtkm::FloatDefault IntersectionHeight = 0.0;
          vtkm::FloatDefault IntersectionWidth = 0.0;

          // data rectangle
          minX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinX);
          maxX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxX);
          minY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinY);
          maxY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxY);

          // caculating intersection of two rectangle (overlapping region)
          minX_intersection = std::max(minX_user, minX_dataset);
          minY_intersection = std::max(minY_user, minY_dataset);
          maxX_intersection = std::min(maxX_user, maxX_dataset);
          maxY_intersection = std::min(maxY_user, maxY_dataset);

          IntersectionHeight = maxY_intersection - minY_intersection;
          IntersectionWidth = maxX_intersection - minX_intersection;

          vtkm::FloatDefault DataArea = (maxX_dataset - minX_dataset) * (maxY_dataset - minY_dataset);

          if ((IntersectionHeight > 0) && (IntersectionWidth > 0) && (minX_intersection < maxX_intersection) && (minY_intersection < maxY_intersection))
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
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputBottomLeft;
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputTopRight;
      };

      class FiberMean : public vtkm::worklet::WorkletMapField
      {
      public:
        FiberMean(const vtkm::Pair<vtkm::Float64, vtkm::Float64> &minAxis,
                  const vtkm::Pair<vtkm::Float64, vtkm::Float64> &maxAxis)
            : InputBottomLeft(minAxis), InputTopRight(maxAxis){};

        // Input and Output Parameters
        using ControlSignature = void(FieldIn, FieldIn, FieldIn, FieldIn, FieldOut);

        using ExecutionSignature = void(_1, _2, _3, _4, _5);
        using InputDomain = _1;

        template <typename MinX,
                  typename MaxX,
                  typename MinY,
                  typename MaxY,
                  typename OutCellFieldType>

        VTKM_EXEC void operator()(const MinX &EnsembleMinX,
                                  const MaxX &EnsembleMaxX,
                                  const MinY &EnsembleMinY,
                                  const MaxY &EnsembleMaxY,
                                  OutCellFieldType &probability) const
        {

          // User defined rectangle(trait)
          vtkm::FloatDefault minX_user = 0.0;
          minX_user = static_cast<vtkm::FloatDefault>(InputBottomLeft.first);
          vtkm::FloatDefault minY_user = 0.0;
          minY_user = static_cast<vtkm::FloatDefault>(InputBottomLeft.second);
          vtkm::FloatDefault maxX_user = 0.0;
          maxX_user = static_cast<vtkm::FloatDefault>(InputTopRight.first);
          vtkm::FloatDefault maxY_user = 0.0;
          maxY_user = static_cast<vtkm::FloatDefault>(InputTopRight.second);

          vtkm::FloatDefault minX_dataset = 0.0;
          vtkm::FloatDefault minY_dataset = 0.0;
          vtkm::FloatDefault maxX_dataset = 0.0;
          vtkm::FloatDefault maxY_dataset = 0.0;

          // data rectangle
          minX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinX);
          maxX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxX);
          minY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinY);
          maxY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxY);

    
          vtkm::FloatDefault Xmean = 0.0;
          vtkm::FloatDefault Ymean = 0.0;
          Xmean = (minX_dataset + maxX_dataset) / 2;
          Ymean = (minY_dataset + maxY_dataset) / 2;

          if ((Xmean <= maxX_user) && (Xmean >= minX_user) && (Ymean <= maxY_user) && (Ymean >= minY_user))
          {
            probability = 1.0;
            return;
          }
          else
          {
            probability = 0.0;
            return;
          }
        }

      private:
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputBottomLeft;
        vtkm::Pair<vtkm::Float64, vtkm::Float64> InputTopRight;
      };

    }
  }
}
#endif
