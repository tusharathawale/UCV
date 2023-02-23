#include "vtkSubsampleUncertaintyEnsemble.h"

#include "SubsampleUncertaintyEnsemble.h"

#include "vtkAOSDataArrayTemplate.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkmlib/ImageDataConverter.h"

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Error.h>

#include <algorithm>
#include <array>

namespace
{

inline std::array<int, 3> ExtentToDimensions(const std::array<int, 6>& wholeExtent)
{
  return { wholeExtent[1] - wholeExtent[0] + 1,
           wholeExtent[3] - wholeExtent[2] + 1,
           wholeExtent[5] - wholeExtent[4] + 1 };
}

inline std::array<int, 6> DimensionsToExtent(
    const std::array<int, 3>& dimensions, const std::array<int, 3>& origin = { 0, 0, 0 })
{
  return { origin[0], dimensions[0] + origin[0] - 1,
           origin[1], dimensions[1] + origin[1] - 1,
           origin[2], dimensions[2] + origin[2] - 1 };
}

template <typename T>
inline vtkSmartPointer<vtkAOSDataArrayTemplate<T>> VtkmToVtkArray(
  const vtkm::cont::ArrayHandle<vtkm::Vec<T, 64>>& vtkmArray)
{
  vtkNew<vtkAOSDataArrayTemplate<T>> vtkArray;
  vtkArray->SetNumberOfComponents(64);

  vtkmArray.SyncControlArray();
  const vtkm::Id size = vtkmArray.GetNumberOfValues() * 64;
  auto bufferInfo = vtkmArray.GetBuffers()[0].GetHostBufferInfo();
  vtkm::cont::internal::TransferredBuffer transfer = bufferInfo.TransferOwnership();
  vtkArray->SetVoidArray(reinterpret_cast<T*>(transfer.Memory),
                         size,
                         0,
                         vtkAbstractArray::VTK_DATA_ARRAY_USER_DEFINED);
  vtkArray->SetArrayFreeFunction(transfer.Delete);

  return vtkArray;
}

} // abstract namespace

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkSubsampleUncertaintyEnsemble);

vtkSubsampleUncertaintyEnsemble::vtkSubsampleUncertaintyEnsemble() = default;

vtkSubsampleUncertaintyEnsemble::~vtkSubsampleUncertaintyEnsemble() = default;

// This filter is going to resize the structured. The VTK streaming pipeline uses
// this information to request parts of the data. Make sure we accurately report
// the size of data that this filter will return.
int vtkSubsampleUncertaintyEnsemble::RequestInformation(vtkInformation* request,
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  this->Superclass::RequestInformation(request, inputVector, outputVector);

  // Get the output info object. The superclass has already copied information from input
  // to output, so we are going to edit the whole extent to the new extent this filter
  // will create.
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  std::array<int, 6> wholeExtent;
  outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent.data());

  std::array<int, 3> dimensions = ExtentToDimensions(wholeExtent);

  for (auto& dim : dimensions)
  {
    dim = (dim + this->BlockSize - 1) / this->BlockSize;
  }

  wholeExtent = DimensionsToExtent(dimensions);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent.data(), 6);

  return 1;
}

// Because this filter resized the grid, it altered the whole extent of the data (see
// RequestInformation). Thus, we need to alter the requested extent to match the whole
// extent of the input.
int vtkSubsampleUncertaintyEnsemble::RequestUpdateExtent(vtkInformation* request,
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  this->Superclass::RequestUpdateExtent(request, inputVector, outputVector);

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

  std::array<int, 6> wholeExtent;
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent.data());
  std::array<int, 3> wholeDimensions = ExtentToDimensions(wholeExtent);

  std::array<int, 6> updateExtent;
  inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), updateExtent.data());
  std::array<int, 3> updateDimensions = ExtentToDimensions(updateExtent);

  for (std::size_t i = 0; i < 3; ++i)
  {
    updateDimensions[i] = std::min(updateDimensions[i] * this->BlockSize, wholeDimensions[i]);
  }

  updateExtent = DimensionsToExtent(updateDimensions,
                                    { wholeExtent[0], wholeExtent[2], wholeExtent[4] });
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), updateExtent.data(), 6);

  return 1;
}

int vtkSubsampleUncertaintyEnsemble::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData* output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  std::array<int, 6> extent;
  input->GetExtent(extent.data());
  std::array<int, 3> dimensions = ExtentToDimensions(extent);
  if (((dimensions[0] % this->BlockSize) != 0) ||
      ((dimensions[1] % this->BlockSize) != 0) ||
      ((dimensions[2] % this->BlockSize) != 0))
  {
    vtkErrorMacro(<< "Dimensions of image must be divisible by the block size.");
    return 0;
  }

  try
  {
    // Convert the input dataset to VTK-m
    vtkm::cont::DataSet in = tovtkm::Convert(input, tovtkm::FieldsFlag::PointsAndCells);

    vtkm::filter::uncertainty::SubsampleUncertaintyEnsemble filter;
    filter.SetMeanSuffix(this->MeanSuffix);
    filter.SetEnsembleSuffix(this->EnsembleSuffix);
    filter.SetBlockSize(this->BlockSize);
    vtkm::cont::DataSet result = filter.Execute(in);

    if (!fromvtkm::Convert(result, output, input))
    {
      vtkErrorMacro(<< "Unable to convert VTK-m DataSet back to VTK.");
      return 0;
    }

    // As of this writing, the VTK-m to VTK convert will not handle arrays with tuples of
    // size 64. This should be resolved once VTK updates to VTK-m 2.1, but until then
    // here is a hack to copy them.
    vtkm::IdComponent numFields = result.GetNumberOfFields();
    for (vtkm::IdComponent fieldIdx = 0; fieldIdx < numFields; ++fieldIdx)
    {
      const vtkm::cont::Field& field = result.GetField(fieldIdx);
      if (field.IsPointField() && !output->GetPointData()->HasArray(field.GetName().c_str()))
      {
        vtkm::cont::UnknownArrayHandle fieldData = field.GetData();
        if (fieldData.IsValueType<vtkm::Vec<float, 64>>())
        {
          auto vtkArray = VtkmToVtkArray(
                fieldData.AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::Vec<float, 64>>>());
          vtkArray->SetName(field.GetName().c_str());
          output->GetPointData()->AddArray(vtkArray);
        }
        else if (fieldData.IsValueType<vtkm::Vec<double, 64>>())
        {
          auto vtkArray = VtkmToVtkArray(
                fieldData.AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::Vec<double, 64>>>());
          vtkArray->SetName(field.GetName().c_str());
          output->GetPointData()->AddArray(vtkArray);
        }
        else
        {
          // Just ignore.
        }
      }
    }
  }
  catch (const vtkm::cont::Error& e)
  {
    vtkErrorMacro(<< "VTK-m error: " << e.GetMessage());
    return 0;
  }

  return 1;
}

void vtkSubsampleUncertaintyEnsemble::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "MeanSuffix: " << this->MeanSuffix << "\n";
  os << indent << "EnsembleSuffix: " << this->EnsembleSuffix << "\n";
  os << indent << "BlockSize: " << this->BlockSize << "\n";
}

VTK_ABI_NAMESPACE_END
