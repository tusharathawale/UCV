#include "vtkUncertainContour.h"

#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include "vtkmlib/DataSetConverters.h"
#include "../../ucvworklet/ExtractingMeanStdev.hpp"
#include "../../ucvworklet/EntropyIndependentGaussian.hpp"

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkUncertainContour);

vtkUncertainContour::vtkUncertainContour() = default;
vtkUncertainContour::~vtkUncertainContour() = default;

std::string vtkUncertainContour::GetInputArrayName(
    int index, vtkInformationVector **inputVector)
{
    int association = this->GetInputArrayAssociation(index, inputVector);
    vtkDataArray *inputArray = this->GetInputArrayToProcess(index, inputVector);
    if ((association != vtkDataObject::FIELD_ASSOCIATION_POINTS) || (inputArray == nullptr))
    {
        vtkErrorMacro("Invalid minimum array; array missing or not a point array.");
        return 0;
    }

    const char *scalarFieldName = inputArray->GetName();
    if (!scalarFieldName || scalarFieldName[0] == '\0')
    {
        scalarFieldName = tovtkm::NoNameVTKFieldName();
    }

    return scalarFieldName;
}

vtkm::cont::DataSet vtkUncertainContour::CallUncertainContourWorklet(vtkm::cont::DataSet vtkmDataSet)
{
    std::cout << "debug CallUncertainContourWorklet" << std::endl;
    std::cout << "---debug isovalue is " << this->IsoValue << std::endl;
    auto outputDataSet = vtkmDataSet;

    // Processing current ensemble data sets based on uncertianty countour
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdevArray;
    auto resolveType = [&](auto &concreteArray)
    {
        vtkm::cont::Invoker invoke;
        vtkm::Id numPoints = concreteArray.GetNumberOfValues();
        std::cout << "---debug numPoints is " << numPoints << std::endl;

        auto concreteArrayView = vtkm::cont::make_ArrayHandleView(concreteArray, 0, numPoints);

        invoke(ExtractingMeanStdevEnsembles{}, concreteArrayView, meanArray, stdevArray);
        printSummary_ArrayHandle(meanArray, std::cout);
        printSummary_ArrayHandle(stdevArray, std::cout);

        vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
        vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
        vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

        invoke(EntropyIndependentGaussian<4, 16>{this->IsoValue}, vtkmDataSet.GetCellSet(), meanArray, stdevArray, crossProbability, numNonZeroProb, entropy);

        outputDataSet.AddCellField(this->ContourProbabilityName, crossProbability);
        outputDataSet.AddCellField(this->NumberNonzeroProbabilityName, numNonZeroProb);
        outputDataSet.AddCellField(this->EntropyName, entropy);

        // vtkm::io::VTKDataSetWriter writeCross(outputFileName);
        // writeCross.WriteDataSet(outputDataSet);
    };

    vtkmDataSet.GetField("ensemble")
        .GetData()
        .CastAndCallWithExtractedArray(resolveType);

    return outputDataSet;
}

int vtkUncertainContour::RequestData(
    vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    // Get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    try
    {
        // Convert the input dataset to VTK-m
        
        std::cout << "original vtk data set" << std::endl;

        std::cout << *input << std::endl;

        // write out data for checking

        vtkm::cont::DataSet in = tovtkm::Convert(input, tovtkm::FieldsFlag::PointsAndCells);

        std::cout << "debug input vtkm data" << std::endl;
        in.PrintSummary(std::cout);

        vtkm::cont::DataSet result = this->CallUncertainContourWorklet(in);

        // Convert the result back from VTKm to VTK.
        // It would be easier if there was a simple method to just convert from general
        // vtkm::cont::DataSet to vtkDataSet. However, that does not exist so you have
        // to copy the vtkDataSet structure in VTK and copy the new fields over. I think
        // it is done this way to prevent creating new arrays for what should be shallow
        // copies. (Maybe in the future the data sharing will be good enough where that
        // is not an issue.)

        output->ShallowCopy(input);
        auto copyField = [&](const std::string &fieldName)
        {
            vtkDataArray *resultingArray = fromvtkm::Convert(result.GetCellField(fieldName));
            if (resultingArray == nullptr)
            {
                vtkWarningMacro(<< "Unable to convert result array " << fieldName << " from VTK-m to VTK.");
                return;
            }
            output->GetCellData()->AddArray(resultingArray);
            resultingArray->FastDelete();
        };
        // Copy the field from the vtkm to vtk and adding associated array into it
        copyField(this->ContourProbabilityName);
        copyField(this->NumberNonzeroProbabilityName);
        copyField(this->EntropyName);
        output->GetCellData()->SetActiveScalars(this->ContourProbabilityName.c_str());
    }
    catch (const vtkm::cont::Error &e)
    {
        vtkErrorMacro(<< "VTK-m error: " << e.GetMessage());
        return 0;
    }

    return 1;
}

void vtkUncertainContour::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);

    os << indent << "IsoValue: " << this->IsoValue << "\n";
    os << indent << "Distribution: " << this->Distribution << "\n";
    os << indent << "ContourProbabilityName: " << this->ContourProbabilityName << "\n";
    os << indent << "NumberNonzeroProbabilityName: " << this->NumberNonzeroProbabilityName << "\n";
    os << indent << "EntropyName: " << this->EntropyName << "\n";
}

VTK_ABI_NAMESPACE_END