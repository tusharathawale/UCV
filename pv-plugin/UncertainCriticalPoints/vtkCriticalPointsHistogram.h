/**
 * @class vtkCriticalPointsHistogram
 * @brief Convert a field representing an ensemble of values into a histogram representation.
 *
 * This filter uses a histogram to model the uncertainty of a field. That is, every value in
 * the field as an associated histogram modeling the probability of the value at that point
 * in the field. This histogram is specified using two field arrays. The first field array
 * has `Vec`s with the number of components equal to the number of bins. The second field array
 * contains the histogram "edges", which are the field values between the bins. The size of
 * the `Vec`s in this array must be one more than the number of histogram bins.
 *
 * The algorithm for this filter works by creating a 3D histogram of values to determine the
 * probable location of critical points. This 3D histogram is independent of the histograms
 * determining the field probabilities, and its size is specified separately.
 */

#ifndef vtkCriticalPointsHistogram_h
#define vtkCriticalPointsHistogram_h

#include "vtkDataObject.h"
#include "vtkDataSetAlgorithm.h"
#include "vtkUncertainCriticalPointsFiltersModule.h" // for export macro
#include "vtkmlib/vtkmInitializer.h"

VTK_ABI_NAMESPACE_BEGIN

class VTKUNCERTAINCRITICALPOINTSFILTERS_EXPORT vtkCriticalPointsHistogram
  : public vtkDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkCriticalPointsHistogram, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkCriticalPointsHistogram* New();

  ///@{
  /// @brief Specify the number of bins to build each histogram to.
  ///
  /// The algorithm of this filter uses a 3D histogram to determine the probable locations of
  /// critical points. The number of computations grows cubically with respect to the number
  /// of bins selected. A larger number of bins provides a more accurate result at the cost
  /// of cubically growing time. The default value is 5.
  vtkSetMacro(NumberOfBins, int);
  vtkGetMacro(NumberOfBins, int);
  ///@}

  /// @brief Specify the name of the field to use for the histogram bin density values.
  ///
  /// The probability of the input field is model using a histogram. The histograms are
  /// represented with two field arrays. The first field array contains, for each value, a
  /// tuple with the fraction of values in the bin (or, equivalently, the probability the
  /// value is in that bin).
  void SetHistogramDensity(const char* name)
  {
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
  }

  /// @brief Specify the name of the field to use for the histogram bin edge values.
  ///
  /// The probability of the input field is model using a histogram. The histograms are
  /// represented with two field arrays. The second field array contains, for each value, a
  /// tuple with the "edge" values for the histogram. The number of edges is one more than
  /// the number of bins. The first edge is the minimum sample value and the last edge is
  /// the maximum edge value. All in between specify the value represented by where two bins
  /// are adjacent.
  void SetHistogramEdges(const char* name)
  {
    this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
  }

  ///@{
  /// @brief Specify the name of the histogram density output field array name.
  vtkSetMacro(MinimumProbabilityName, std::string);
  vtkGetMacro(MinimumProbabilityName, std::string);
  ///@}

protected:
  vtkCriticalPointsHistogram();
  ~vtkCriticalPointsHistogram();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  int NumberOfBins = 5;
  std::string MinimumProbabilityName = "MinimumProbability";

private:
  vtkCriticalPointsHistogram(const vtkCriticalPointsHistogram&) = delete;
  void operator=(const vtkCriticalPointsHistogram&) = delete;

  vtkmInitializer Initializer;

  std::string GetInputArrayName(int index, vtkInformationVector** inputVector);
};

VTK_ABI_NAMESPACE_END

#endif // vtkCriticalPointsHistogram_h
