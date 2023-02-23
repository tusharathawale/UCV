#include "vtkContourUncertainEnsemble.h"

#include "ContourUncertainEnsemble.h"

#include "vtkAOSDataArrayTemplate.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include "vtkmlib/DataArrayConverters.h"
#include "vtkmlib/DataSetConverters.h"

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkContourUncertainEnsemble);

vtkContourUncertainEnsemble::vtkContourUncertainEnsemble() = default;
vtkContourUncertainEnsemble::~vtkContourUncertainEnsemble() = default;

std::string vtkContourUncertainEnsemble::GetInputArrayName(
  int index, vtkInformationVector** inputVector)
{
  int association = this->GetInputArrayAssociation(index, inputVector);
  vtkDataArray* inputArray = this->GetInputArrayToProcess(index, inputVector);
  if ((association != vtkDataObject::FIELD_ASSOCIATION_POINTS) || (inputArray == nullptr))
  {
    vtkErrorMacro("Invalid minimum array; array missing or not a point array.");
    return 0;
  }

  const char* scalarFieldName = inputArray->GetName();
  if (!scalarFieldName || scalarFieldName[0] == '\0')
  {
    scalarFieldName = tovtkm::NoNameVTKFieldName();
  }

  return scalarFieldName;
}

void vtkContourUncertainEnsemble::SetInputMean(const char* name)
{
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
}
void vtkContourUncertainEnsemble::SetInputMean(int fieldAttributeType)
{
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, fieldAttributeType);
}
void vtkContourUncertainEnsemble::SetInputEnsemble(const char* name)
{
  this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, name);
}
void vtkContourUncertainEnsemble::SetInputEnsemble(int fieldAttributeType)
{
  this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, fieldAttributeType);
}

int vtkContourUncertainEnsemble::RequestData(
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

    // We are expecting the ensemble field to have tuple sizes of 64. This is currently
    // not directly supported and results in errors. This should be handled correctly
    // once we move to VTK-m 2.1, but for now hack the movement.
    vtkDataArray* ensembleArray = this->GetInputArrayToProcess(1, inputVector);
    if (ensembleArray->GetNumberOfComponents() == 64)
    {
      bool converted = false;
      auto tryVtkToVtkm = [&](auto vtkArray) {
        if (vtkArray == nullptr)
        {
          return;
        }
        auto vtkmArray =
          tovtkm::DataArrayToArrayHandle<std::remove_pointer_t<decltype(vtkArray)>, 64>::Wrap(vtkArray);
        in.AddPointField(vtkArray->GetName(), vtkmArray);
        converted = true;
      };
      tryVtkToVtkm(vtkAOSDataArrayTemplate<float>::SafeDownCast(ensembleArray));
      tryVtkToVtkm(vtkAOSDataArrayTemplate<double>::SafeDownCast(ensembleArray));
      if (!converted)
      {
        vtkErrorMacro(<< "Unable to convert ensemble field.");
        return 0;
      }
    }

    vtkm::filter::uncertainty::ContourUncertainEnsemble filter;
    filter.SetIsoValue(this->IsoValue);
    filter.SetCrossProbabilityName(this->ContourProbabilityName);
    filter.SetNumberNonzeroProbabilityName(this->NumberNonzeroProbabilityName);
    filter.SetEntropyName(this->EntropyName);
    filter.SetMeanField(this->GetInputArrayName(0, inputVector));
    filter.SetEnsembleField(this->GetInputArrayName(1, inputVector));
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
      vtkDataArray* resultingArray = fromvtkm::Convert(result.GetCellField(fieldName));
      if (resultingArray == nullptr)
      {
        vtkWarningMacro(<< "Unable to convert result array " << fieldName << " from VTK-m to VTK.");
        return;
      }
      output->GetCellData()->AddArray(resultingArray);
      resultingArray->FastDelete();
    };
    copyField(this->ContourProbabilityName);
    copyField(this->NumberNonzeroProbabilityName);
    copyField(this->EntropyName);
    output->GetCellData()->SetActiveScalars(this->ContourProbabilityName.c_str());
  }
  catch (const vtkm::cont::Error& e)
  {
    vtkErrorMacro(<< "VTK-m error: " << e.GetMessage());
    return 0;
  }

  return 1;
}

void vtkContourUncertainEnsemble::PrintSelf(ostream&os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "IsoValue: " << this->IsoValue << "\n";
  os << indent << "ContourProbabilityName: " << this->ContourProbabilityName << "\n";
  os << indent << "NumberNonzeroProbabilityName: " << this->NumberNonzeroProbabilityName << "\n";
  os << indent << "EntropyName: " << this->EntropyName << "\n";
}

VTK_ABI_NAMESPACE_END
