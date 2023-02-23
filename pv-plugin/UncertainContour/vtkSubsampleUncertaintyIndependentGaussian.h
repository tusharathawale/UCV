/**
 * @class vtkSubsampleUncertaintyIndependentGaussian
 * @brief Subsamples a regular grid and captures the uncertainty as a normal distribution.
 *
 * This filter takes a uniform grid and subsamples it by a particular factor.
 * As the grid is subsampled, the uncertainty is captured by considering the
 * block of points that get reduced to the single point. This uncertainty is
 * modeled as a normal distirbution with a single independent Gaussian function.
 * It is captured as the mean and standard deviation within the block.
 */

#ifndef vtkSubsampleUncertaintyIndependentGaussian_h
#define vtkSubsampleUncertaintyIndependentGaussian_h

#include "vtkUncertainContourFiltersModule.h" // for export macro
#include "vtkImageAlgorithm.h"
#include "vtkmlib/vtkmInitializer.h"

VTK_ABI_NAMESPACE_BEGIN
class VTKUNCERTAINCONTOURFILTERS_EXPORT vtkSubsampleUncertaintyIndependentGaussian
  : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkSubsampleUncertaintyIndependentGaussian, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkSubsampleUncertaintyIndependentGaussian *New();

  ///@{
  /// \brief The suffix used for fields modeling uncertainty.
  ///
  /// When this filter subsamples data, it models the uncertainty of each field by finding
  /// the mean and standard deviation in each region that is grouped together. These are captured
  /// as a pair of fields in the output. They are given the name of the input field appended
  /// with the appropriate suffix.
  ///
  vtkSetMacro(MeanSuffix, std::string);
  vtkGetMacro(MeanSuffix, std::string);
  vtkSetMacro(StdevSuffix, std::string);
  vtkGetMacro(StdevSuffix, std::string);
  ///@}

  ///@{
  /// \brief Specifies the reduction factor for the subsampling
  ///
  vtkSetMacro(BlockSize, int);
  vtkGetMacro(BlockSize, int);
  ///@}

protected:
  vtkSubsampleUncertaintyIndependentGaussian();
  ~vtkSubsampleUncertaintyIndependentGaussian() override;

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  std::string MeanSuffix = "_mean";
  std::string StdevSuffix = "_stdev";
  int BlockSize = 4;

private:
  vtkSubsampleUncertaintyIndependentGaussian(
    const vtkSubsampleUncertaintyIndependentGaussian&) = delete;
  void operator=(const vtkSubsampleUncertaintyIndependentGaussian&) = delete;

  vtkmInitializer Initializer;
};
VTK_ABI_NAMESPACE_END

#endif //vtkSubsampleUncertaintyIndependentGaussian_h
