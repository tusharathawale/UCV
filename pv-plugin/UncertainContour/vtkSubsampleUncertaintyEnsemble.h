/**
 * @class vtkSubsampleUncertaintyEnsemble
 * @brief Subsamples a regular grid and captures the uncertainty as an ensemble.
 *
 * This filter takes a uniform grid and subsamples it by a particular factor.
 * As the grid is subsampled, the uncertainty is captured by considering the
 * block of points that get reduced to the single point. This uncertainty is
 * modeled as the mean value and an ensemble collection, which captures all the
 * values for the block.
 */

#ifndef vtkSubsampleUncertaintyEnsemble_h
#define vtkSubsampleUncertaintyEnsemble_h

#include "vtkUncertainContourFiltersModule.h" // for export macro
#include "vtkImageAlgorithm.h"
#include "vtkmlib/vtkmInitializer.h"

VTK_ABI_NAMESPACE_BEGIN
class VTKUNCERTAINCONTOURFILTERS_EXPORT vtkSubsampleUncertaintyEnsemble
  : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkSubsampleUncertaintyEnsemble, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkSubsampleUncertaintyEnsemble *New();

  ///@{
  /// \brief The suffix used for fields modeling uncertainty.
  ///
  /// When this filter subsamples data, it models the uncertainty of each field by finding
  /// the mean and an ensemble of the original values. These are captured as a pair of fields
  /// in the output. They are given the name of the input field appended with the appropriate
  /// suffix.
  ///
  vtkSetMacro(MeanSuffix, std::string);
  vtkGetMacro(MeanSuffix, std::string);
  vtkSetMacro(EnsembleSuffix, std::string);
  vtkGetMacro(EnsembleSuffix, std::string);
  ///@}

  ///@{
  /// \brief Specifies the reduction factor for the subsampling
  ///
  vtkSetMacro(BlockSize, int);
  vtkGetMacro(BlockSize, int);
  ///@}

protected:
  vtkSubsampleUncertaintyEnsemble();
  ~vtkSubsampleUncertaintyEnsemble() override;

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  std::string MeanSuffix = "_mean";
  std::string EnsembleSuffix = "_ensemble";
  int BlockSize = 4;

private:
  vtkSubsampleUncertaintyEnsemble(const vtkSubsampleUncertaintyEnsemble&) = delete;
  void operator=(const vtkSubsampleUncertaintyEnsemble&) = delete;

  vtkmInitializer Initializer;
};
VTK_ABI_NAMESPACE_END

#endif //vtkSubsampleUncertaintyEnsemble_h
