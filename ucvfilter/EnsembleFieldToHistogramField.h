#ifndef vtk_m_filter_uncertainty_EnsembleFieldToHistogramField
#define vtk_m_filter_uncertainty_EnsembleFieldToHistogramField

#include <vtkm/filter/FilterField.h>

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

/// @brief Convert a field representing an ensemble of values into a histogram representation.
///
/// One way to estimate the uncertainty of a system is to collect an ensemble of values. For
/// example, one might run a simulation multiple times with varying initial conditions. This
/// filter takes a group of ensemble measurements that have been placed into a single field.
/// The field contains `Vec` values with one component for each member of the ensemble. The
/// filter converts this representation to a histogram representation of the uncertainty. The
/// histogram provides an estimate of the probability distribution function (PDF) for the
/// members of the uncertainty. The value of each histogram bin is the probability that an
/// ensemble member is located in that bin.
///
/// The histograms for the field values are represented by two field arrays. The first array
/// is the histogram "density," which is the probability values for the bins. This array has
/// `Vec`s the size of _N_ where _N_ is the number of bins selected for the histogram. The
/// second array contains the histogram "edges," which define the scalar values that each
/// bin represents. This array has `Vec`s the size of _N_+1. The first value is the minimum
/// scalar value; the last value is the maximum value, and the remaining values provide the
/// subsequent scalar values in between histogram bins. They are evenly spaced with
/// (max-min)/_N_ distance between them.
///
/// Set the active field to specify which field to convert from an ensemble representation
/// to a histogram representation.
class EnsembleFieldToHistogramField : public vtkm::filter::FilterField
{
  vtkm::IdComponent NumberOfBins = 0;
  std::string HistogramEdgesName = "HistogramEdges";

public:
  VTKM_CONT EnsembleFieldToHistogramField();

  /// @brief Specify the number of bins to build each histogram to.
  ///
  /// A larger number of bins provides a closer fit to the ensemble data. However, larger
  /// numbers require more memory to represent. Also, using too many bins will result in
  /// an "overfitting" that assumes the ensemble represents all samples that could possibly
  /// happen.
  ///
  /// If this value is unset or set to less than 2, Sturge's rule will be used to select
  /// a number of bins based on the size of the input ensemble.
  VTKM_CONT void SetNumberOfBins(vtkm::IdComponent numBins) { this->NumberOfBins = numBins; }
  /// @copydoc SetNumberOfBins
  VTKM_CONT vtkm::IdComponent GetNumberOfBins() const { return this->NumberOfBins; }

  /// @brief Specify the name of the histogram density output field array name.
  VTKM_CONT void SetHistogramDensityName(const std::string& name)
  {
    this->SetOutputFieldName(name);
  }
  /// @copydoc SetHistogramDensityName
  VTKM_CONT std::string GetHistogramDensityName() const { return this->GetOutputFieldName(); }

  /// @brief Specify the name of the histogram edges output field array name.
  VTKM_CONT void SetHistogramEdgesName(const std::string& name) { this->HistogramEdgesName = name; }
  /// @copydoc SetHistogramEdgesName
  VTKM_CONT std::string GetHistogramEdgesName() const { return this->HistogramEdgesName; }

protected:
  VTKM_CONT vtkm::cont::DataSet DoExecute(const vtkm::cont::DataSet& inData) override;
};

}
}
} // namespace vtkm::filter::uncertainty

#endif //vtk_m_filter_uncertainty_EnsembleFieldToHistogramField
