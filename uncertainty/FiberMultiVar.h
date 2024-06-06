//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_uncertainty_FiberMultiVar_h
#define vtk_m_filter_uncertainty_FiberMultiVar_h

#include <vtkm/filter/FilterField.h>
#include <vtkm/filter/uncertainty/vtkm_filter_uncertainty_export.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{
class VTKM_FILTER_UNCERTAINTY_EXPORT FiberMultiVar : public vtkm::filter::FilterField
{
  vtkm::Vec<vtkm::Float64, 3> bottomLeft;
  vtkm::Vec<vtkm::Float64, 3> topRight;

public:

  // x1, y1, z1 
  VTKM_CONT void SetBottomLeftAxis(const vtkm::Vec<vtkm::Float64, 3>& bottomLeftCoordinate)
  {
    this->bottomLeft = bottomLeftCoordinate;
  }
  
  // x2, y2, z2 
  VTKM_CONT void SetTopRightAxis(const vtkm::Vec<vtkm::Float64, 3>& topRightCoordinate)
  {
    this->topRight = topRightCoordinate;
  }

  // (x1, x2)
  VTKM_CONT void SetMinX(const std::string& fieldName)
  {
    this->SetActiveField(0, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxX(const std::string& fieldName)
  {
    this->SetActiveField(1, fieldName, vtkm::cont::Field::Association::Points);
  }
  
  // (y1, y2)
  VTKM_CONT void SetMinY(const std::string& fieldName)
  {
    this->SetActiveField(2, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxY(const std::string& fieldName)
  {
    this->SetActiveField(3, fieldName, vtkm::cont::Field::Association::Points);
  }

  // (z1, z2)
  VTKM_CONT void SetMinZ(const std::string& fieldName)
  {
    this->SetActiveField(4, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxZ(const std::string& fieldName)
  {
    this->SetActiveField(5, fieldName, vtkm::cont::Field::Association::Points);
  }

private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& input) override;
};
}
}
}
#endif
