//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_filter_uncertainty_ContourUncertainEnsemble2D_h
#define vtk_m_filter_uncertainty_ContourUncertainEnsemble2D_h

#include <vtkm/filter/FilterField.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

/// \brief Computes the probability of the location of a contour for 2d.
///
/// This filter computes the probable location of a contour of a field
/// with uncertainty. Unlike a typical contour filter that extracts a polygonal surface,
/// this filter simply writes a cell field giving the probability of a contour being in
/// that cell.
///
/// This filter takes as input a field with a `Vec` containing the ensemble members
/// it use the mulivariate gaussian to compute the uncertainty isocountour
/// this filter is supposed to be merged together with ContourUncertainEnsemble in future
/// it compute the mean in the worklet instead of using a explicit mean array as input parameter
class ContourUncertainEnsemble2D : public vtkm::filter::FilterField
{
  std::string NumberNonzeroProbabilityName = "num_nonzero_probability";
  std::string EntropyName = "entropy";
  vtkm::Float64 IsoValue = 0.1;

public:
  VTKM_CONT ContourUncertainEnsemble2D();
  VTKM_CONT ~ContourUncertainEnsemble2D() = default;
  ///@{
  /// Specifies the fields capturing ensemble values and mean of the data.
  VTKM_CONT void SetEnsembleField(const std::string& fieldName)
  {
    this->SetActiveField(0, fieldName, vtkm::cont::Field::Association::Points);
  }
  ///@}

  ///@{
  /// Specifies the contour value.
  VTKM_CONT void SetIsoValue(vtkm::Float64 value) { this->IsoValue = value; }
  VTKM_CONT vtkm::Float64 GetIsoValue() const { return this->IsoValue; }
  ///@}

  ///@{
  /// Specifies the name of the output field that captures the probability of the contour existing
  /// in each cell.
  VTKM_CONT void SetCrossProbabilityName(const std::string& name)
  {
    this->SetOutputFieldName(name);
  }
  VTKM_CONT const std::string& GetCrossProbabilityName() const
  {
    return this->GetOutputFieldName();
  }
  ///@}

  ///@{
  /// Specifies the name of the output field that captures the number of possible marching
  /// contour cases for each cell.
  ///
  /// When running a marching contour algorithm on certain data, an exact contour case
  /// can be determined. But when marching uncertain data, multiple cases are possible.
  /// This field reports the number of cases that have a non-zero probability. The more
  /// cases that are possible, the more uncertain the contour.
  ///
  VTKM_CONT void SetNumberNonzeroProbabilityName(const std::string& name)
  {
    this->NumberNonzeroProbabilityName = name;
  }
  VTKM_CONT const std::string& GetNumberNonzeroProbabilityName() const
  {
    return this->NumberNonzeroProbabilityName;
  }
  ///@}

  ///@{
  /// Specifies the name of the output field that captures the entropy of the possible
  /// marching contour cases for each cell.
  ///
  /// When running a marching contour algorithm on certain data, an exact contour case
  /// can be determined. But when marching uncertain data, multiple cases are possible.
  /// This field reports the entropy of the contour cases. A high entropy means that
  /// the case random, and thus the actual contour is uncertain.
  ///
  VTKM_CONT void SetEntropyName(const std::string& name) { this->EntropyName = name; }
  VTKM_CONT const std::string& GetEntropyName() const { return this->EntropyName; }
  ///@}

protected:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& input) override;
};

}
}
} // namespace vtkm::filter::uncertainty

#endif //vtk_m_filter_uncertainty_ContourUncertainEnsemble2D_h
