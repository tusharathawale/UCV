/**
 * @class vtkSubsampleUncertaintyUniform
 * @brief Subsamples a regular grid and captures the uncertainty as a uniform distribution.
 *
 * This filter takes a uniform grid and subsamples it by a particular factor.
 * As the grid is subsampled, the uncertainty is captured by considering the
 * block of points that get reduced to the single point. This uncertainty is
 * modeled as a uniform distribution and captured as the minimum and maximum
 * value within the block.
 */

#ifndef vtkSubsampleUncertaintyUniform_h
#define vtkSubsampleUncertaintyUniform_h

#include "vtkUncertainContourFiltersModule.h" // for export macro
#include "vtkImageAlgorithm.h"
#include "vtkmlib/vtkmInitializer.h"

VTK_ABI_NAMESPACE_BEGIN
class VTKUNCERTAINCONTOURFILTERS_EXPORT vtkSubsampleUncertaintyUniform : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkSubsampleUncertaintyUniform, vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkSubsampleUncertaintyUniform *New();

  ///@{
  /// \brief The suffix used for fields modeling uncertainty.
  ///
  /// When this filter subsamples data, it models the uncertainty of each field by finding
  /// the minimum and maximum value in each region that is grouped together. These are captured
  /// as a pair of fields in the output. They are given the name of the input field appended
  /// with the appropriate suffix.
  ///
  vtkSetMacro(MinSuffix, std::string);
  vtkGetMacro(MinSuffix, std::string);
  vtkSetMacro(MaxSuffix, std::string);
  vtkGetMacro(MaxSuffix, std::string);
  ///@}

  ///@{
  /// \brief Specifies the reduction factor for the subsampling
  ///
  vtkSetMacro(BlockSize, int);
  vtkGetMacro(BlockSize, int);
  ///@}

protected:
  vtkSubsampleUncertaintyUniform();
  ~vtkSubsampleUncertaintyUniform() override;

  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestUpdateExtent(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  std::string MinSuffix = "_min";
  std::string MaxSuffix = "_max";
  int BlockSize = 4;

private:
  vtkSubsampleUncertaintyUniform(const vtkSubsampleUncertaintyUniform&) = delete;
  void operator=(const vtkSubsampleUncertaintyUniform&) = delete;

  vtkmInitializer Initializer;
};
VTK_ABI_NAMESPACE_END

#endif //vtkSubsampleUncertaintyUniform_h
