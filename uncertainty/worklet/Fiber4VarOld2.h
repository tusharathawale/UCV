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
#include <vtkm/worklet/WorkletPointNeighborhood.h>

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
      class FiberMonteCarlo : public vtkm::worklet::WorkletPointNeighborhood
      {
      public:
        // Worklet Input
        // Fiber(const std::vector<std::pair<double, double>>& minAxis,
        //      const std::vector<std::pair<double, double>>& maxAxis)
        //  : SetMinAxis(minAxis), SetMaxAxis(maxAxis){};
        FiberMonteCarlo(const vtkm::Vec<vtkm::Float64, 4> &minAxis,
                        const vtkm::Vec<vtkm::Float64, 4> &maxAxis,
                        const vtkm::Id numSamples)
            : SetMinAxis(minAxis), SetMaxAxis(maxAxis), NumSamples(numSamples){};

        // Input and Output Parameters
        using ControlSignature = void(CellSetIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn,  FieldOut);
        using ExecutionSignature = void(_2, _3, _4, _5, _6, _7, _8, _9, _10);
        using InputDomain = _1;

        // Template
        template <typename MinX,
                  typename MaxX,
                  typename MinY,
                  typename MaxY,
                  typename MinZ,
                  typename MaxZ,
                  typename MinW,
                  typename MaxW,
                  typename OutCellFieldType>
        // Operator
        VTKM_EXEC void operator()(const MinX &EnsembleMinX,
                                  const MaxX &EnsembleMaxX,
                                  const MinY &EnsembleMinY,
                                  const MaxY &EnsembleMaxY,
                                  const MinZ &EnsembleMinZ,
                                  const MaxZ &EnsembleMaxZ,
                                  const MinW &EnsembleMinW,
                                  const MaxW &EnsembleMaxW,
                                  OutCellFieldType &probability) const
        {
          // User defined rectangle(trait)
          vtkm::FloatDefault minX_user = 0.0;
          minX_user = static_cast<vtkm::FloatDefault>(SetMinAxis[0]);
          vtkm::FloatDefault minY_user = 0.0;
          minY_user = static_cast<vtkm::FloatDefault>(SetMinAxis[1]);
          vtkm::FloatDefault minZ_user = 0.0;
          minZ_user = static_cast<vtkm::FloatDefault>(SetMinAxis[2]);
          vtkm::FloatDefault minW_user = 0.0;
          minW_user = static_cast<vtkm::FloatDefault>(SetMinAxis[3]);
          vtkm::FloatDefault maxX_user = 0.0;
          maxX_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[0]);
          vtkm::FloatDefault maxY_user = 0.0;
          maxY_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[1]);
          vtkm::FloatDefault maxZ_user = 0.0;
          maxZ_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[2]);
          vtkm::FloatDefault maxW_user = 0.0;
          maxW_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[3]);
          // vtkm::FloatDefault TraitArea = (maxX_user - minX_user) * (maxY_user - minY_user);

          vtkm::FloatDefault minX_dataset = 0.0;
          vtkm::FloatDefault minY_dataset = 0.0;
          vtkm::FloatDefault maxX_dataset = 0.0;
          vtkm::FloatDefault maxY_dataset = 0.0;
          vtkm::FloatDefault minZ_dataset = 0.0;
          vtkm::FloatDefault maxZ_dataset = 0.0;
          vtkm::FloatDefault minW_dataset = 0.0;
          vtkm::FloatDefault maxW_dataset = 0.0;

          vtkm::FloatDefault N1 = 0.0;
          vtkm::FloatDefault N2 = 0.0;
          vtkm::FloatDefault N3 = 0.0;
          vtkm::FloatDefault N4 = 0.0;
          vtkm::IdComponent NonZeroCases = 0;
          vtkm::FloatDefault MCProbability = 0.0;

          // data rectangle
          minX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinX);
          maxX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxX);
          minY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinY);
          maxY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxY);
          minZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinZ);
          maxZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxZ);
          minW_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinW);
          maxW_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxW);


          // if data rectangle is zero, there is no uncertainty, return zero
          // TODO needs 4var implementation
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
          std::uniform_real_distribution<vtkm::FloatDefault> GenerateN3(minZ_dataset, maxZ_dataset);
          std::uniform_real_distribution<vtkm::FloatDefault> GenerateN4(minW_dataset, maxW_dataset);

          for (vtkm::IdComponent i = 0; i < this->NumSamples; i++)
          {
            N1 = GenerateN1(gen);
            N2 = GenerateN2(gen);
            N3 = GenerateN3(gen);
            N4 = GenerateN4(gen);
            // increase the case number when the data is located in user rectangle
            if ((N1 > minX_user) && (N1 < maxX_user) && (N2 > minY_user) && (N2 < maxY_user) && (N3 > minZ_user) && (N3 < maxZ_user) && (N4 > minW_user) && (N4 < maxW_user))
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
        vtkm::Vec<vtkm::Float64, 4> SetMinAxis;
        vtkm::Vec<vtkm::Float64, 4> SetMaxAxis;
        vtkm::Id NumSamples = 1;
      };

      class FiberClosedForm : public vtkm::worklet::WorkletPointNeighborhood
      {
      public:
        // Worklet Input
        // Fiber(const std::vector<std::pair<double, double>>& minAxis,
        //      const std::vector<std::pair<double, double>>& maxAxis)
        //  : SetMinAxis(minAxis), SetMaxAxis(maxAxis){};
        FiberClosedForm(const vtkm::Vec<vtkm::Float64, 4> &minAxis,
                        const vtkm::Vec<vtkm::Float64, 4> &maxAxis)
            : SetMinAxis(minAxis), SetMaxAxis(maxAxis){};

        // Input and Output Parameters
        using ControlSignature = void(CellSetIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn,  FieldOut);
        using ExecutionSignature = void(_2, _3, _4, _5, _6, _7, _8, _9, _10);
        using InputDomain = _1;

        // Template
        template <typename MinX,
                  typename MaxX,
                  typename MinY,
                  typename MaxY,
                  typename MinZ,
                  typename MaxZ,
                  typename MinW,
                  typename MaxW,
                  typename OutCellFieldType>
        // Operator
        VTKM_EXEC void operator()(const MinX &EnsembleMinX,
                                  const MaxX &EnsembleMaxX,
                                  const MinY &EnsembleMinY,
                                  const MaxY &EnsembleMaxY,
                                  const MinZ &EnsembleMinZ,
                                  const MaxZ &EnsembleMaxZ,
                                  const MinW &EnsembleMinW,
                                  const MaxW &EnsembleMaxW,
                                  OutCellFieldType &probability) const
        {
          // User defined rectangle(trait)
          vtkm::FloatDefault minX_user = 0.0;
          minX_user = static_cast<vtkm::FloatDefault>(SetMinAxis[0]);
          vtkm::FloatDefault minY_user = 0.0;
          minY_user = static_cast<vtkm::FloatDefault>(SetMinAxis[1]);
          vtkm::FloatDefault minZ_user = 0.0;
          minZ_user = static_cast<vtkm::FloatDefault>(SetMinAxis[2]);
          vtkm::FloatDefault minW_user = 0.0;
          minW_user = static_cast<vtkm::FloatDefault>(SetMinAxis[3]);
          vtkm::FloatDefault maxX_user = 0.0;
          maxX_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[0]);
          vtkm::FloatDefault maxY_user = 0.0;
          maxY_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[1]);
          vtkm::FloatDefault maxZ_user = 0.0;
          maxZ_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[2]);
          vtkm::FloatDefault maxW_user = 0.0;
          maxW_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[3]);


          // data rectangle
          vtkm::FloatDefault minX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinX);
          vtkm::FloatDefault maxX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxX);
          vtkm::FloatDefault minY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinY);
          vtkm::FloatDefault maxY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxY);
          vtkm::FloatDefault minZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinZ);
          vtkm::FloatDefault maxZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxZ);
          vtkm::FloatDefault minW_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinW);
          vtkm::FloatDefault maxW_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxW);

          // caculating intersection area
          vtkm::FloatDefault minX_intersection = std::max(minX_user, minX_dataset);
          vtkm::FloatDefault minY_intersection = std::max(minY_user, minY_dataset);
          vtkm::FloatDefault maxX_intersection = std::min(maxX_user, maxX_dataset);
          vtkm::FloatDefault maxY_intersection = std::min(maxY_user, maxY_dataset);
          vtkm::FloatDefault minZ_intersection = std::max(minZ_user, minZ_dataset);
          vtkm::FloatDefault maxZ_intersection = std::min(maxZ_user, maxZ_dataset);
          vtkm::FloatDefault minW_intersection = std::max(minW_user, minW_dataset);
          vtkm::FloatDefault maxW_intersection = std::min(maxW_user, maxW_dataset);



          if ((minX_intersection < maxX_intersection) and (minY_intersection < maxY_intersection) and (minZ_intersection < maxZ_intersection) and (minW_intersection < maxW_intersection))
          {
            vtkm::FloatDefault volume_intersection = (maxX_intersection - minX_intersection) * (maxY_intersection - minY_intersection) * (maxZ_intersection - minZ_intersection) * (maxW_intersection - minW_intersection);
            vtkm::FloatDefault volume_data = (maxX_dataset - minX_dataset) * (maxY_dataset - minY_dataset) * (maxZ_dataset - minZ_dataset) * (maxW_dataset - minW_dataset);
            probability = volume_intersection / volume_data;
          }
          else
          {
            probability = 0.0;
          }
        }

      private:
        vtkm::Vec<vtkm::Float64, 4> SetMinAxis;
        vtkm::Vec<vtkm::Float64, 4> SetMaxAxis;
      };

      class FiberMean : public vtkm::worklet::WorkletPointNeighborhood
      {
      public:
        FiberMean(const vtkm::Vec<vtkm::Float64, 4> &minAxis,
                  const vtkm::Vec<vtkm::Float64, 4> &maxAxis)
            : SetMinAxis(minAxis), SetMaxAxis(maxAxis){};

        // Input and Output Parameters
        using ControlSignature = void(CellSetIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn,  FieldOut);
        using ExecutionSignature = void(_2, _3, _4, _5, _6, _7, _8, _9, _10);
        using InputDomain = _1;

        // Template
        template <typename MinX,
                  typename MaxX,
                  typename MinY,
                  typename MaxY,
                  typename MinZ,
                  typename MaxZ,
                  typename MinW,
                  typename MaxW,
                  typename OutCellFieldType>
        // Operator
        VTKM_EXEC void operator()(const MinX &EnsembleMinX,
                                  const MaxX &EnsembleMaxX,
                                  const MinY &EnsembleMinY,
                                  const MaxY &EnsembleMaxY,
                                  const MinZ &EnsembleMinZ,
                                  const MaxZ &EnsembleMaxZ,
                                  const MinW &EnsembleMinW,
                                  const MaxW &EnsembleMaxW,
                                  OutCellFieldType &probability) const
        {
          // User defined rectangle(trait)
          vtkm::FloatDefault minX_user = 0.0;
          minX_user = static_cast<vtkm::FloatDefault>(SetMinAxis[0]);
          vtkm::FloatDefault minY_user = 0.0;
          minY_user = static_cast<vtkm::FloatDefault>(SetMinAxis[1]);
          vtkm::FloatDefault minZ_user = 0.0;
          minZ_user = static_cast<vtkm::FloatDefault>(SetMinAxis[2]);
          vtkm::FloatDefault minW_user = 0.0;
          minW_user = static_cast<vtkm::FloatDefault>(SetMinAxis[3]);
          vtkm::FloatDefault maxX_user = 0.0;
          maxX_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[0]);
          vtkm::FloatDefault maxY_user = 0.0;
          maxY_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[1]);
          vtkm::FloatDefault maxZ_user = 0.0;
          maxZ_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[2]);
          vtkm::FloatDefault maxW_user = 0.0;
          maxW_user = static_cast<vtkm::FloatDefault>(SetMaxAxis[3]);


          // data rectangle
          vtkm::FloatDefault minX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinX);
          vtkm::FloatDefault maxX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxX);
          vtkm::FloatDefault minY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinY);
          vtkm::FloatDefault maxY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxY);
          vtkm::FloatDefault minZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinZ);
          vtkm::FloatDefault maxZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxZ);
          vtkm::FloatDefault minW_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinW);
          vtkm::FloatDefault maxW_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxW);


          vtkm::FloatDefault Xmean = 0.0;
          vtkm::FloatDefault Ymean = 0.0;
          vtkm::FloatDefault Zmean = 0.0;
          vtkm::FloatDefault Wmean = 0.0;
          Xmean = (minX_dataset + maxX_dataset) / 2;
          Ymean = (minY_dataset + maxY_dataset) / 2;
          Zmean = (minZ_dataset + maxZ_dataset) / 2;
          Wmean = (minW_dataset + maxW_dataset) / 2;

          if ((Xmean <= maxX_user) && (Xmean >= minX_user) && (Ymean <= maxY_user) && (Ymean >= minY_user) && (Zmean <= maxZ_user) && (Zmean >= minZ_user) && (Wmean <= maxW_user) && (Wmean >= minW_user))
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
        vtkm::Vec<vtkm::Float64, 4> SetMinAxis;
        vtkm::Vec<vtkm::Float64, 4> SetMaxAxis;
      };

    }
  }
}
#endif
