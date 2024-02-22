//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_uncertainty_Fiber_h
#define vtk_m_filter_uncertainty_Fiber_h

#include <vtkm/filter/FilterField.h>
#include <vtkm/filter/uncertainty/vtkm_filter_uncertainty_export.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{
class VTKM_FILTER_UNCERTAINTY_EXPORT Fiber : public vtkm::filter::FilterField
{
  vtkm::Pair<vtkm::Float64, vtkm::Float64> minAxis;
  vtkm::Pair<vtkm::Float64, vtkm::Float64> maxAxis;

public:
  VTKM_CONT void SetMinAxis(const vtkm::Pair<vtkm::Float64, vtkm::Float64>& minCoordinate)
  {
    this->minAxis = minCoordinate;
  }

  VTKM_CONT void SetMaxAxis(const vtkm::Pair<vtkm::Float64, vtkm::Float64>& maxCoordinate)
  {
    this->maxAxis = maxCoordinate;
  }
  VTKM_CONT void SetMinOne(const std::string& fieldName)
  {
    this->SetActiveField(0, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxOne(const std::string& fieldName)
  {
    this->SetActiveField(1, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMinTwo(const std::string& fieldName)
  {
    this->SetActiveField(2, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxTwo(const std::string& fieldName)
  {
    this->SetActiveField(3, fieldName, vtkm::cont::Field::Association::Points);
  }

private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& input) override;
};
}
}
}
#endif
