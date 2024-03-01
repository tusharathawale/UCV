#include "vtkFiberSurfaceUncertain.h"

#include "../../Fiber.h"

#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"

#include "vtkmlib/DataSetConverters.h"

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkFiberSurfaceUncertain);

vtkFiberSurfaceUncertain::vtkFiberSurfaceUncertain() = default;
vtkFiberSurfaceUncertain::~vtkFiberSurfaceUncertain() = default;

std::string vtkFiberSurfaceUncertain::GetInputArrayName(
  int index, vtkInformationVector** inputVector)
{
  int association = this->GetInputArrayAssociation(index, inputVector);
  vtkDataArray* inputArray = this->GetInputArrayToProcess(index, inputVector);
  if ((association != vtkDataObject::FIELD_ASSOCIATION_POINTS) || (inputArray == nullptr))
  {
    vtkErrorMacro("Invalid array; array missing or not a point array.");
    return 0;
  }

  const char* scalarFieldName = inputArray->GetName();
  if (!scalarFieldName || scalarFieldName[0] == '\0')
  {
    scalarFieldName = tovtkm::NoNameVTKFieldName();
  }

  return scalarFieldName;
}

int vtkFiberSurfaceUncertain::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkDataSet* input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet* output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  try
  {
    // Convert the input dataset to VTK-m
    vtkm::cont::DataSet in = tovtkm::Convert(input, tovtkm::FieldsFlag::PointsAndCells);

    vtkm::filter::uncertainty::Fiber filter;

    filter.SetMinAxis({ this->Axis1Minimum, this->Axis2Minimum} );
    filter.SetMaxAxis({ this->Axis1Maximum, this->Axis2Maximum} );

    filter.SetMinOne(this->GetInputArrayName(0, inputVector));
    filter.SetMaxOne(this->GetInputArrayName(1, inputVector));
    filter.SetMinTwo(this->GetInputArrayName(2, inputVector));
    filter.SetMaxTwo(this->GetInputArrayName(3, inputVector));

    filter.SetNumSamples(this->NumSamples);

    switch (this->Approach)
    {
      case MONTE_CARLO: filter.SetApproach("MonteCarlo"); break;
      case CLOSED_FORM: filter.SetApproach("ClosedForm"); break;
      default:
        vtkErrorMacro("Invalid Approach.");
        return 0;
    }

    vtkm::cont::DataSet result = filter.Execute(in);

    // Convert the result back.
    // It would be easier if there was a simple method to just convert from general
    // vtkm::cont::DataSet to vtkDataSet. However, that does not exist so you have
    // to copy the vtkDataSet structure in VTK and copy the new fields over. I think
    // it is done this way to prevent creating new arrays for what should be shallow
    // copies. (Maybe in the future the data sharing will be good enough where that
    // is not an issue.)
    output->ShallowCopy(input);
    auto copyField = [&](const std::string& fieldName) {
      vtkDataArray* resultingArray = fromvtkm::Convert(result.GetPointField(fieldName));
      if (resultingArray == nullptr)
      {
        vtkWarningMacro(<< "Unable to convert result array " << fieldName << " from VTK-m to VTK.");
        return;
      }
      output->GetPointData()->AddArray(resultingArray);
      resultingArray->FastDelete();
    };
    copyField("OutputProbability");
    output->GetPointData()->SetActiveScalars("OutputProbability");
  }
  catch (const vtkm::cont::Error& e)
  {
    vtkErrorMacro(<< "VTK-m error: " << e.GetMessage());
    return 0;
  }

  return 1;
}

void vtkFiberSurfaceUncertain::PrintSelf(ostream&os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Axis1Minimum: " << this->Axis1Minimum << "\n";
  os << indent << "Axis1Maximum: " << this->Axis1Maximum << "\n";
  os << indent << "Axis2Minimum: " << this->Axis2Minimum << "\n";
  os << indent << "Axis2Maximum: " << this->Axis2Maximum << "\n";
  os << index << "Approach: ";
  switch (this->Approach)
  {
    case MONTE_CARLO: os << "Monte Carlo\n"; break;
    case CLOSED_FORM: os << "Closed Form\n"; break;
    default: os << "?\n"; break;
  }
  os << index << "NumSamples: " << this->NumSamples << "\n";
}

VTK_ABI_NAMESPACE_END
