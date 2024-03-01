/**
 * @class vtkFiberSurfaceUncertain
 * @brief Finds the probability of fiber surface location given a uniform distribution of uncertainty.
 *
 * This filter computes the probable location of a fiber surface with uncertainty.
 * Unlike a typical fiber surface filter that extracts a polygonal surface, this
 * filter simply writes a point field giving the probability of a surface
 * being inside the closed fiber surface.
 *
 * This filter uses a uniform distribution model for the uncertainty of a
 * field. To use the filter, you must specify a field with the minimum value
 * and one with the maximum value.
 */

#ifndef vtkFiberSurfaceUncertain_h
#define vtkFiberSurfaceUncertain_h

#include "vtkUncertainFiberSurfaceFiltersModule.h" // for export macro
#include "vtkDataObject.h"
#include "vtkDataSetAlgorithm.h"
#include "vtkmlib/vtkmInitializer.h"

VTK_ABI_NAMESPACE_BEGIN

class VTKUNCERTAINFIBERSURFACEFILTERS_EXPORT vtkFiberSurfaceUncertain : public vtkDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkFiberSurfaceUncertain, vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  static vtkFiberSurfaceUncertain *New();

  ///@{
  /// @brief The input fields used to model uncertainty.
  ///
  /// When this filter determines the probable locations of a contour,
  /// it uses a uniform distribution model of the uncertainty. This model
  /// is defined by a minimum value and maximum value. These values are
  /// expected to be stored in separate fields, which are selected here.
  ///
  /// Note that these are convenience methods that call `SetInputArrayToProcess`
  /// to set input 0 for the minimum and input 1 for the maximum.
  ///
  void SetInput1Minimum(const char *name)
  {
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
  }
  void SetInput1Minimum(int fieldAttributeType)
  {
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, fieldAttributeType);
  }
  void SetInput1Maximum(const char *name)
  {
    this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
  }
  void SetInput1Maximum(int fieldAttributeType)
  {
    this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, fieldAttributeType);
  }
  void SetInput2Minimum(const char *name)
  {
    this->SetInputArrayToProcess(2, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
  }
  void SetInput2Minimum(int fieldAttributeType)
  {
    this->SetInputArrayToProcess(2, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, fieldAttributeType);
  }
  void SetInput2Maximum(const char *name)
  {
    this->SetInputArrayToProcess(3, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
  }
  void SetInput2Maximum(int fieldAttributeType)
  {
    this->SetInputArrayToProcess(3, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, fieldAttributeType);
  }
  ///@}

  ///@{
  /// @brief Specify the bounds of the fiber surface in field space.
  ///
  /// A fiber surface in general draws a closed curve in the 2D space representing
  /// all pairs formed from 2 fields. The fiber surface contains all points containing
  /// field pairs on this curve. Currently, this filter only supports an axis aligned
  /// bounding box for the curve, which is defined by the minimum and maximum value
  /// for each field.
  vtkSetMacro(Axis1Minimum, double);
  vtkGetMacro(Axis1Minimum, double);
  vtkSetMacro(Axis1Maximum, double);
  vtkGetMacro(Axis1Maximum, double);
  vtkSetMacro(Axis2Minimum, double);
  vtkGetMacro(Axis2Minimum, double);
  vtkSetMacro(Axis2Maximum, double);
  vtkGetMacro(Axis2Maximum, double);
  ///@}

  enum {
    MONTE_CARLO = 0,
    CLOSED_FORM = 1
  };

  ///@{
  /// @brief Specify the approach used to compute the probabilistic fiber surface.
  vtkSetMacro(Approach, int);
  vtkGetMacro(Approach, int);
  ///@}

  ///@{
  /// @brief Specify the number of samples used for Monte Carlo computation.
  vtkSetMacro(NumSamples, int);
  vtkGetMacro(NumSamples, int);
  ///@}

protected:
  vtkFiberSurfaceUncertain();
  ~vtkFiberSurfaceUncertain();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  double Axis1Minimum = 0.0;
  double Axis1Maximum = 0.0;
  double Axis2Minimum = 0.0;
  double Axis2Maximum = 0.0;

  int Approach = MONTE_CARLO;

  int NumSamples = 500;

private:
  vtkFiberSurfaceUncertain(const vtkFiberSurfaceUncertain&) = delete;
  void operator=(const vtkFiberSurfaceUncertain&) = delete;

  vtkmInitializer Initializer;

  std::string GetInputArrayName(int index, vtkInformationVector** inputVector);
};

VTK_ABI_NAMESPACE_END

#endif //vtkFiberSurfaceUncertain_h
