//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_filter_uncertainty_Fiber_separate_h
#define vtk_m_filter_uncertainty_Fiber_separate_h

#include <vtkm/filter/Filter.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{
class Fiber3Var : public vtkm::filter::Filter
{
  vtkm::Vec<vtkm::Float64, 3> minAxis;
  vtkm::Vec<vtkm::Float64, 3> maxAxis;
  std::string Approach = "ClosedForm"; //MonteCarlo or ClosedForm
  vtkm::Id NumSamples = 500;

public:
  VTKM_CONT void SetMinAxis(const vtkm::Vec<vtkm::Float64, 3>& minCoordinate)
  {
    this->minAxis = minCoordinate;
  }

  VTKM_CONT void SetMaxAxis(const vtkm::Vec<vtkm::Float64, 3>& maxCoordinate)
  {
    this->maxAxis = maxCoordinate;
  }
  VTKM_CONT void SetMinX(const std::string& fieldName)
  {
    this->SetActiveField(0, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxX(const std::string& fieldName)
  {
    this->SetActiveField(1, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMinY(const std::string& fieldName)
  {
    this->SetActiveField(2, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxY(const std::string& fieldName)
  {
    this->SetActiveField(3, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMinZ(const std::string& fieldName)
  {
    this->SetActiveField(4, fieldName, vtkm::cont::Field::Association::Points);
  }
  VTKM_CONT void SetMaxZ(const std::string& fieldName)
  {
    this->SetActiveField(5, fieldName, vtkm::cont::Field::Association::Points);
  }

  VTKM_CONT void SetNumSamples(const vtkm::Id& numSamples)
  {
    this->NumSamples=numSamples;
  }

   VTKM_CONT void SetApproach(const std::string& approach)
  {
    this->Approach=approach;
  } 


private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& input) override;
};
}
}
}
#endif
