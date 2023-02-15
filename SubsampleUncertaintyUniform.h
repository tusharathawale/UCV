//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_filter_uncertainty_SubsampleUncertaintyUniform_h
#define vtk_m_filter_uncertainty_SubsampleUncertaintyUniform_h

#include <vtkm/filter/Filter.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

/// \brief Subsamples a regular grid and captures the uncertainty as a uniform distribution.
/// This uncertainty is modeled using the minimum and maximum value.
class SubsampleUncertaintyUniform : public vtkm::filter::Filter
{
  std::string MinSuffix = "_min";
  std::string MaxSuffix = "_max";
  vtkm::IdComponent BlockSize = 4;

public:
  SubsampleUncertaintyUniform() = default;

  ///@{
  /// \brief The suffix used for fields modeling uncertainty.
  ///
  /// When this filter subsamples data, it models the uncertainty of each field by finding
  /// the minimum and maximum value in each region that is grouped together. These are captured
  /// as a pair of fields in the output. They are given the name of the input field appended
  /// with the appropriate suffix.
  ///
  VTKM_CONT void SetMinSuffix(const std::string& minSuffix) { this->MinSuffix = minSuffix; }
  VTKM_CONT const std::string& GetMinSuffix() const { return this->MinSuffix; }
  VTKM_CONT void SetMaxSuffix(const std::string& maxSuffix) { this->MaxSuffix = maxSuffix; }
  VTKM_CONT const std::string& GetMaxSuffix() const { return this->MaxSuffix; }
  ///@}

  ///@{
  /// \brief Specifies the reduction factor for the subsampling
  ///
  VTKM_CONT void SetBlockSize(vtkm::IdComponent blocksize) { this->BlockSize = blocksize; }
  VTKM_CONT vtkm::IdComponent GetBlockSize() const { return this->BlockSize; }
  ///@}

private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& input) override;
};

}
}
} // namespace vtkm::filter::uncertainty

#endif //vtk_m_filter_uncertainty_SubsampleUncertaintyUniform_h
