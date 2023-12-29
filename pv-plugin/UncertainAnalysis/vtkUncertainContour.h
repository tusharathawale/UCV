/**
 * @class vtkUncertainContour
 * @brief Finds the probability of contour location given data with uncertainty.
 *
 * This filter computes the probable location of a field with uncertainty.
 * Unlike a typical contour filter that extracts a polygonal surface, this
 * filter simply writes a cell field giving the probability of a contour
 * being in that cell.
 *
 */

#ifndef vtkUncertainContour_h
#define vtkUncertainContour_h

#include "vtkUncertainAnalysisFiltersModule.h" // for export macro
//#include "vtkImageAlgorithm.h" // Currently only support regular grids. Superclass may change.
#include "vtkDataSetAlgorithm.h"
#include "vtkmlib/vtkmInitializer.h"
#include "vtkmlib/DataSetConverters.h"

VTK_ABI_NAMESPACE_BEGIN
class VTKUNCERTAINANALYSISFILTERS_EXPORT vtkUncertainContour : public vtkDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkUncertainContour, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkUncertainContour *New();

  ///@{
  /// \brief The scalar value to use for the isosurface.
  ///
  /// As with any contouring algorithm, a scalar value is specified to extract
  /// the surface where the field is equal to this value. (Of course, this
  /// filter does not actually extract the surface but rather determines
  /// the probability of the surface being in each cell.)
  ///
  vtkSetMacro(IsoValue, float);
  vtkGetMacro(IsoValue, float);
  ///@}

  ///@{
  /// \brief The distribution used for computing the crossing probability.
  ///
  ///
  vtkSetMacro(Distribution, std::string);
  vtkGetMacro(Distribution, std::string);
  ///@}

  ///@{
  /// \brief The num of samples used for img
  ///
  ///
  vtkSetMacro(NumSamples, int);
  vtkGetMacro(NumSamples, int);
  ///@}

  ///@{
  /// The name of the output field giving the probability of the contour existing in each cell.
  ///
  vtkSetMacro(ContourProbabilityName, std::string);
  vtkGetMacro(ContourProbabilityName, std::string);
  ///@}

  ///@{
  /// The name of the output field giving the number of possible marching
  /// contour cases for each cell.
  ///
  /// When running a marching contour algorithm on certain data, an exact contour case
  /// can be determined. But when marching uncertain data, multiple cases are possible.
  /// This field reports the number of cases that have a non-zero probability. The more
  /// cases that are possible, the more uncertain the contour.
  ///
  vtkSetMacro(NumberNonzeroProbabilityName, std::string);
  vtkGetMacro(NumberNonzeroProbabilityName, std::string);
  ///@}

  ///@{
  /// The name of the output field giving the entropy of the possible
  /// marching contour cases for each cell.
  ///
  /// When running a marching contour algorithm on certain data, an exact contour case
  /// can be determined. But when marching uncertain data, multiple cases are possible.
  /// This field reports the entropy of the contour cases. A high entropy means that
  /// the case random, and thus the actual contour is uncertain.
  ///
  vtkSetMacro(EntropyName, std::string);
  vtkGetMacro(EntropyName, std::string);
  ///@}

  vtkm::cont::DataSet CallUncertainContourWorklet(vtkm::cont::DataSet input);

protected:
  vtkUncertainContour();
  ~vtkUncertainContour();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  double IsoValue = 0.0;
  std::string ContourProbabilityName = "contour_probability";
  std::string NumberNonzeroProbabilityName = "num_nonzero_probability";
  std::string EntropyName = "entropy";
  std::string Distribution = "mvg"; // ig (indepedent gaussian) uni (uniform distribution) kde (kernel density estimation)
  int NumSamples = 1000;

private:
  vtkUncertainContour(const vtkUncertainContour&) = delete;
  void operator=(const vtkUncertainContour&) = delete;

  vtkmInitializer Initializer;

  std::string GetInputArrayName(int index, vtkInformationVector** inputVector);
};

VTK_ABI_NAMESPACE_END

#endif //vtkUncertainContour_h