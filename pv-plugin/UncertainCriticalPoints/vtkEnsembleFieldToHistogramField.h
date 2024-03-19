/**
 * @class vtkEnsembleFieldToHistogramField
 * @brief Convert a field representing an ensemble of values into a histogram representation.
 *
 * One way to estimate the uncertainty of a system is to collect an ensemble of values. For
 * example, one might run a simulation multiple times with varying initial conditions. This
 * filter takes a group of ensemble measurements that have been placed into a single field.
 * The field contains `Vec` values with one component for each member of the ensemble. The
 * filter converts this representation to a histogram representation of the uncertainty. The
 * histogram provides an estimate of the probability distribution function (PDF) for the
 * members of the uncertainty. The value of each histogram bin is the probability that an
 * ensemble member is located in that bin.
 *
 * The histograms for the field values are represented by two field arrays. The first array
 * is the histogram "density," which is the probability values for the bins. This array has
 * `Vec`s the size of _N_ where _N_ is the number of bins selected for the histogram. The
 * second array contains the histogram "edges," which define the scalar values that each
 * bin represents. This array has `Vec`s the size of _N_+1. The first value is the minimum
 * scalar value; the last value is the maximum value, and the remaining values provide the
 * subsequent scalar values in between histogram bins. They are evenly spaced with
 * (max-min)/_N_ distance between them.
 *
 * Set the active field to specify which field to convert from an ensemble representation
 * to a histogram representation.
 */

#ifndef vtkEnsembleFieldToHistogramField_h
#define vtkEnsembleFieldToHistogramField_h

#include "vtkDataObject.h"
#include "vtkDataSetAlgorithm.h"
#include "vtkUncertainCriticalPointsFiltersModule.h" // for export macro
#include "vtkmlib/vtkmInitializer.h"

VTK_ABI_NAMESPACE_BEGIN

class VTKUNCERTAINCRITICALPOINTSFILTERS_EXPORT vtkEnsembleFieldToHistogramField
  : public vtkDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkEnsembleFieldToHistogramField, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkEnsembleFieldToHistogramField* New();

  ///@{
  /// @brief Specify the number of bins to build each histogram to.
  ///
  /// A larger number of bins provides a closer fit to the ensemble data. However, larger
  /// numbers require more memory to represent. Also, using too many bins will result in
  /// an "overfitting" that assumes the ensemble represents all samples that could possibly
  /// happen.
  ///
  /// If this value is unset or set to less than 2, Sturge's rule will be used to select
  /// a number of bins based on the size of the input ensemble.
  vtkSetMacro(NumberOfBins, int);
  vtkGetMacro(NumberOfBins, int);
  ///@}

  ///@{
  /// @brief Specify the name of the histogram density output field array name.
  vtkSetMacro(HistogramDensityName, std::string);
  vtkGetMacro(HistogramDensityName, std::string);
  ///@}

  ///@{
  /// @brief Specify the name of the histogram edges output field array name.
  vtkSetMacro(HistogramEdgesName, std::string);
  vtkGetMacro(HistogramEdgesName, std::string);
  ///@}

protected:
  vtkEnsembleFieldToHistogramField();
  ~vtkEnsembleFieldToHistogramField();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  int NumberOfBins = 0;
  std::string HistogramDensityName = "HistogramDensity";
  std::string HistogramEdgesName = "HistogramEdges";

private:
  vtkEnsembleFieldToHistogramField(const vtkEnsembleFieldToHistogramField&) = delete;
  void operator=(const vtkEnsembleFieldToHistogramField&) = delete;

  vtkmInitializer Initializer;

  std::string GetInputArrayName(int index, vtkInformationVector** inputVector);
};

VTK_ABI_NAMESPACE_END

#endif // vtkEnsembleFieldToHistogramField_h
