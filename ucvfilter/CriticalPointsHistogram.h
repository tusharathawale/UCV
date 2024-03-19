#ifndef vtk_m_filter_uncertainty_CriticalPointsHistogram_h
#define vtk_m_filter_uncertainty_CriticalPointsHistogram_h

#include <vtkm/filter/FilterField.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

/// @brief Compute the probable locations of critical points in an uncertain field.
///
/// This filter uses a histogram to model the uncertainty of a field. That is, every value in
/// the field as an associated histogram modeling the probability of the value at that point
/// in the field. This histogram is specified using two field arrays. The first field array
/// has `Vec`s with the number of components equal to the number of bins. The second field array
/// contains the histogram "edges", which are the field values between the bins. The size of
/// the `Vec`s in this array must be one more than the number of histogram bins.
///
/// The algorithm for this filter works by creating a 3D histogram of values to determine the
/// probable location of critical points. This 3D histogram is independent of the histograms
/// determining the field probabilities, and its size is specified separately.
class CriticalPointsHistogram : public vtkm::filter::FilterField
{
  vtkm::IdComponent NumberOfBins = 5;

public:
  VTKM_CONT CriticalPointsHistogram();

  /// @brief Specify the number of bins to use when computing the probabilitic critical points.
  ///
  /// The algorithm of this filter uses a 3D histogram to determine the probable locations of
  /// critical points. The number of computations grows cubically with respect to the number
  /// of bins selected. A larger number of bins provides a more accurate result at the cost
  /// of cubically growing time. The default value is 5.
  VTKM_CONT void SetNumberOfBins(vtkm::IdComponent numBins) { this->NumberOfBins = numBins; }
  /// @copydoc SetNumberOfBins
  VTKM_CONT vtkm::IdComponent GetNumberOfBins() const { return this->NumberOfBins; }

  /// @brief Specify the name of the field to use for the histogram bin density values.
  ///
  /// The probability of the input field is model using a histogram. The histograms are
  /// represented with two field arrays. The first field array contains, for each value, a
  /// `Vec` with the fraction of values in the bin (or, equivalently, the probability the
  /// value is in that bin).
  VTKM_CONT void SetHistogramDensityName(const std::string& name)
  {
    this->SetActiveField(0, name, vtkm::cont::Field::Association::Points);
  }
  /// @copydoc SetHistogramDensityName
  VTKM_CONT std::string GetHistogramDensityName() const { return this->GetActiveFieldName(0); }

  /// @brief Specify the name of the field to use for the histogram bin edge values.
  ///
  /// The probability of the input field is model using a histogram. The histograms are
  /// represented with two field arrays. The second field array contains, for each value, a
  /// `Vec` with the "edge" values for the histogram. The number of edges is one more than
  /// the number of bins. The first edge is the minimum sample value and the last edge is
  /// the maximum edge value. All in between specify the value represented by where two bins
  /// are adjacent.
  VTKM_CONT void SetHistogramEdgesName(const std::string& name)
  {
    this->SetActiveField(1, name, vtkm::cont::Field::Association::Points);
  }
  /// @copydoc SetHistogramEdgesName
  VTKM_CONT std::string GetHistogramEdgesName() const { return this->GetActiveFieldName(1); }

  /// @brief Specify the name of the output field holding the probability of a local minimum value.
  VTKM_CONT void SetMinimumProbabilityName(const std::string& fieldName)
  {
    this->SetOutputFieldName(fieldName);
  }
  // @copydoc SetMinimumProbabilityName
  VTKM_CONT std::string GetMinimumProbabilityName() const { return this->GetOutputFieldName(); }

protected:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& inData) override;
};

}
}
} // namespace vtkm::filter::uncertainty

#endif //vtk_m_filter_uncertainty_CriticalPointsHistogram_h
