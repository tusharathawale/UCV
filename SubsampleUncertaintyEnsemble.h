//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================
#ifndef vtk_m_filter_uncertainty_SubsampleUncertaintyEnsemble_h
#define vtk_m_filter_uncertainty_SubsampleUncertaintyEnsemble_h

#include <vtkm/filter/Filter.h>

namespace vtkm
{
namespace worklet
{
// Forward declaration
template <typename T>
class Keys;
}
}

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

/// \brief Subsamples a regular grid and captures the uncertainty as an enseble.
///
/// In the resulting mesh, the fields are adjusted to capture the mean of each block
/// reduced to a single point. Additionally, a field is added that captures all the
/// original data as a `Vec`.
class SubsampleUncertaintyEnsemble : public vtkm::filter::Filter
{
  std::string EnsembleSuffix = "_ensemble";
  vtkm::IdComponent BlockSize = 4;

public:
  SubsampleUncertaintyEnsemble() = default;

  ///@{
  /// \brief The suffix used for fields modeling ensemble.
  ///
  /// When this filter subsamples data, it captures the original data as an ensemble field. This
  /// new field is given the name of the input field appended with the appropriate suffix.
  ///
  VTKM_CONT void SetEnsembleSuffix(const std::string& suffix) { this->EnsembleSuffix = suffix; }
  VTKM_CONT const std::string& GetEnsembleSuffix() const { return this->EnsembleSuffix; }
  ///@}

  ///@{
  /// \brief Specifies the reduction factor for the subsampling
  ///
  VTKM_CONT void SetBlockSize(vtkm::IdComponent blocksize) { this->BlockSize = blocksize; }
  VTKM_CONT vtkm::IdComponent GetBlockSize() const { return this->BlockSize; }
  ///@}

private:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& input) override;

  VTKM_CONT void MapField(vtkm::cont::DataSet& data,
                          const vtkm::cont::Field& field,
                          const vtkm::worklet::Keys<vtkm::Id>& keys) const;
};

}
}
} // namespace vtkm::filter::uncertainty

#endif //vtk_m_filter_uncertainty_SubsampleUncertaintyEnsemble_h
