//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_filter_uncertainty_SubsampleUncertaintyIndependentGaussian_h
#define vtk_m_filter_uncertainty_SubsampleUncertaintyIndependentGaussian_h

#include <vtkm/filter/Filter.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

/// \brief Subsamples a regular grid and captures the uncertainty as an independent gaussian distribution.
/// This uncertainty is modeled using the mean and standard deviation of the data grouped.
class SubsampleUncertaintyIndependentGaussian : public vtkm::filter::Filter
{
  std::string MeanSuffix = "_mean";
  std::string StdevSuffix = "_stdev";
  vtkm::IdComponent BlockSize = 4;

public:
  SubsampleUncertaintyIndependentGaussian() = default;

  ///@{
  /// \brief The suffix used for fields modeling uncertainty.
  ///
  /// When this filter subsamples data, it models the uncertainty of each field by finding
  /// the mean and standard deviation in each region that is grouped together. These are captured
  /// as a pair of fields in the output. They are given the name of the input field appended
  /// with the appropriate suffix.
  ///
  VTKM_CONT void SetMeanSuffix(const std::string& suffix) { this->MeanSuffix = suffix; }
  VTKM_CONT const std::string& GetMeanSuffix() const { return this->MeanSuffix; }
  VTKM_CONT void SetStdevSuffix(const std::string& suffix) { this->StdevSuffix = suffix; }
  VTKM_CONT const std::string& GetStdevSuffix() const { return this->StdevSuffix; }
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

#endif //vtk_m_filter_uncertainty_SubsampleUncertaintyIndependentGaussian_h
