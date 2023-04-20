/**
 * @class vtkContourUncertainEnsemble2D
 * @brief Finds the probability of contour location given a uniform distribution of uncertainty.
 *
 * This filter computes the probable location of a field with uncertainty.
 * Unlike a typical contour filter that extracts a polygonal surface, this
 * filter simply writes a cell field giving the probability of a contour
 * being in that cell.
 *
 * This filter uses an ensemble model for the uncertainty of a
 * field for 2d data, this is the vtk filter, it calles the vktm filter
 */

#ifndef vtkContourUncertainEnsemble2D_h
#define vtkContourUncertainEnsemble2D_h

#include "vtkUncertainContourFiltersModule.h" // for export macro
#include "vtkImageAlgorithm.h" // Currently only support regular grids. Superclass may change.
#include "vtkmlib/vtkmInitializer.h"

VTK_ABI_NAMESPACE_BEGIN

class VTKUNCERTAINCONTOURFILTERS_EXPORT vtkContourUncertainEnsemble2D : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkContourUncertainEnsemble2D, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkContourUncertainEnsemble2D *New();

  ///@{
  /// \brief The input fields used to model uncertainty.
  ///
  /// When this filter determines the probable locations of a contour,
  /// it uses an ensemble model of the uncertainty. This model
  /// is defined by a mean value and all the members of the ensemble. These values are
  /// expected to be stored in separate fields, which are selected here.
  ///
  /// Note that these are convenience methods that call `SetInputArrayToProcess`
  /// to set input 0 for the mean and input 1 for the standard deviation.
  ///
  void SetInputEnsemble(const char* name);
  void SetInputEnsemble(int fieldAttributeType);
  ///@}

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

protected:
  vtkContourUncertainEnsemble2D();
  ~vtkContourUncertainEnsemble2D();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  double IsoValue = 0.0;
  std::string ContourProbabilityName = "contour_probability";
  std::string NumberNonzeroProbabilityName = "num_nonzero_probability";
  std::string EntropyName = "entropy";

private:
  vtkContourUncertainEnsemble2D(const vtkContourUncertainEnsemble2D&) = delete;
  void operator=(const vtkContourUncertainEnsemble2D&) = delete;

  vtkmInitializer Initializer;

  std::string GetInputArrayName(int index, vtkInformationVector** inputVector);
};

VTK_ABI_NAMESPACE_END

#endif //vtkContourUncertainEnsemble2D_h
