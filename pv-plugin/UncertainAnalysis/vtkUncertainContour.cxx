#include "vtkUncertainContour.h"

#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include "vtkmlib/DataSetConverters.h"
#include "../../ucvworklet/ExtractingMeanStdev.hpp"
#include "../../ucvworklet/EntropyIndependentGaussian.hpp"
#include "../../ucvworklet/ComputeDiffSum.hpp"
#include <vtkm/cont/Algorithm.h>
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

vtkm::cont::DataSet vtkUncertainContour::CallUncertainContourWorklet(vtkm::cont::DataSet vtkmDataSet, vtkm::Id &totalNumEnsemble)
{
    std::cout << "debug CallUncertainContourWorklet" << std::endl;
    std::cout << "---debug isovalue is " << this->IsoValue << std::endl;

    // todo, the new data set only need to get the cellset of the original data
    // auto outputDataSet = vtkmDataSet;

    vtkm::cont::DataSet outputDataSet;
    outputDataSet.SetCellSet(vtkmDataSet.GetCellSet());

    // merging ensemble data sets firstly then call functor
    // get number of ens members in data sets, get all fields and check their fields
    std::string ensSuffix = "ensemble_";
    totalNumEnsemble = 0;
    vtkm::IdComponent numFields = vtkmDataSet.GetNumberOfFields();

    // compute the numFields with ensSuffix in it
    vtkm::IdComponent fieldNameWithEns = 0;
    for (vtkm::IdComponent fieldIndex = 0; fieldIndex < numFields; fieldIndex++)
    {
        auto field =vtkmDataSet.GetField(fieldIndex);
        std::string fieldName = field.GetName();
        if (fieldName.find(ensSuffix) != std::string::npos)
        {
            fieldNameWithEns++;
        }
    }
    std::string nameOfFirstEns = "";
    std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> ensFieldArrays;
    //make sure the ensemble is loaded with sequence through 0 to n
    for (vtkm::IdComponent fieldIndex = 0; fieldIndex < fieldNameWithEns; fieldIndex++)
    {
        std::string fieldName = ensSuffix + std::to_string(fieldIndex);
        auto field = vtkmDataSet.GetField(fieldName);
        std::cout << "load field with field name " << fieldName << std::endl;
        if (fieldName.find(ensSuffix) != std::string::npos)
        {
            // this is one of the ensembles
            totalNumEnsemble++;
            if (nameOfFirstEns == "")
            {
                nameOfFirstEns = fieldName;
            }
            vtkm::cont::ArrayHandle<vtkm::Float64> fieldArray;
            field.GetData().AsArrayHandle(fieldArray);
            ensFieldArrays.push_back(fieldArray);
        }
    }

    std::cout << " debug ok to store vtkm array with explicit type" << std::endl;

    // get size of first ens element (how many point in each entry)
    vtkm::Id lengthOfEnsField = vtkmDataSet.GetField(nameOfFirstEns).GetNumberOfValues();

    // for allEnsArray
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::Float64>
        allEnsemblesArray(totalNumEnsemble);
    allEnsemblesArray.Allocate(lengthOfEnsField);

    auto allEnsWritePortal = allEnsemblesArray.WritePortal();
    // cache the read portal
    using ReadPortalType = typename vtkm::cont::ArrayHandle<vtkm::Float64>::ReadPortalType;
    std::vector<ReadPortalType> ReadPortalList;
    for (int ensIndex = 0; ensIndex < totalNumEnsemble; ensIndex++)
    {
        auto readPortal = ensFieldArrays[ensIndex].ReadPortal();
        ReadPortalList.push_back(readPortal);
    }

    for (int fieldValueIndex = 0; fieldValueIndex < lengthOfEnsField; fieldValueIndex++)
    {
        auto vecValue = allEnsWritePortal.Get(fieldValueIndex);
        for (int ensIndex = 0; ensIndex < totalNumEnsemble; ensIndex++)
        {
            vecValue[ensIndex] = ReadPortalList[ensIndex].Get(fieldValueIndex);
        }
    }

    std::cout << "debug, print runtime allEnsemblesArray" << std::endl;
    printSummary_ArrayHandle(allEnsemblesArray, std::cout);

    vtkm::cont::ArrayHandle<vtkm::Float64> meanArray;
    vtkm::cont::ArrayHandle<vtkm::Float64> stdevArray;

    vtkm::cont::Invoker invoke;
    invoke(ExtractingMeanStdevEnsembles{}, allEnsemblesArray, meanArray, stdevArray);
    // printSummary_ArrayHandle(meanArray, std::cout);
    // printSummary_ArrayHandle(stdevArray, std::cout);

    vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

    // TODO, checking the type of cell to determine if call 4, 16 or 3, 8 as input parameter
    vtkm::IdComponent numberOfPointsPerCell = vtkmDataSet.GetCellSet().GetNumberOfPointsInCell(0);
    if (numberOfPointsPerCell == 4)
    {
        invoke(EntropyIndependentGaussian<4, 16>{this->IsoValue}, vtkmDataSet.GetCellSet(), meanArray, stdevArray, crossProbability, numNonZeroProb, entropy);
    }
    else if (numberOfPointsPerCell == 3)
    {
        invoke(EntropyIndependentGaussian<3, 8>{this->IsoValue}, vtkmDataSet.GetCellSet(), meanArray, stdevArray, crossProbability, numNonZeroProb, entropy);
    }
    else
    {
        throw std::runtime_error("numberOfPointsPerCell is " + std::to_string(numberOfPointsPerCell) + " is unsupported");
    }

    outputDataSet.AddCellField(this->ContourProbabilityName, crossProbability);
    outputDataSet.AddCellField(this->NumberNonzeroProbabilityName, numNonZeroProb);
    outputDataSet.AddCellField(this->EntropyName, entropy);

    // For element without ens, create the runtime array without the specific element
    // the code is similar to the place creating allEnsemblesArray
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::Float64>
        withoutOneEnsembleArray(totalNumEnsemble - 1);
    withoutOneEnsembleArray.Allocate(lengthOfEnsField);
    auto withoutOneEnsWritePortal = withoutOneEnsembleArray.WritePortal();
    // cache the read portal
    using ReadPortalType = typename vtkm::cont::ArrayHandle<vtkm::Float64>::ReadPortalType;

    // Results with no ens for each ensemble
    // for each case with out the specific ensemble
    for (int noEnsId = 0; noEnsId < totalNumEnsemble; noEnsId++)
    {
        // init value for no_ens computation
        vtkm::cont::ArrayHandle<vtkm::Float64> noEnscrossProbability;
        vtkm::cont::ArrayHandle<vtkm::Id> noEnsnumNonZeroProb;
        vtkm::cont::ArrayHandle<vtkm::Float64> noEnsentropy;
        vtkm::cont::ArrayHandle<vtkm::Float64> noEnsmeanArray;
        vtkm::cont::ArrayHandle<vtkm::Float64> noEnsstdevArray;

        // the first loop is each value in the lengthOfEnsField
        // the second loop is each ensemble value
        for (vtkm::Id pointIndex = 0; pointIndex < lengthOfEnsField; pointIndex++)
        {
            // the noEnsId should be skipped at each iteration
            int currEnsId = 0;   // currEnsId is from 0 to totalNumEnsemble-2
            int actualEnsId = 0; // actualEnsId is from 0 to totalNumEnsemble-1
            while (actualEnsId < totalNumEnsemble)
            {
                if (actualEnsId == noEnsId)
                {
                    // skip this ens ensemble info
                    actualEnsId++;
                    continue;
                }
                auto vecValue = withoutOneEnsWritePortal.Get(pointIndex);
                vecValue[currEnsId] = ensFieldArrays[actualEnsId].ReadPortal().Get(pointIndex);
                actualEnsId++;
                currEnsId++;
            }
        }

        if (noEnsId == 0)
        {
            std::cout << "debug non ens id 0 vec value " << withoutOneEnsWritePortal.Get(0)[0] << std::endl;
        }

        // compute ens for this case
        // compute mean and stdev
        // std::cout << "process the no ens with id " << noEnsId << std::endl;
        invoke(ExtractingMeanStdevEnsembles{}, withoutOneEnsembleArray, noEnsmeanArray, noEnsstdevArray);
        vtkm::IdComponent numberOfPointsPerCell = vtkmDataSet.GetCellSet().GetNumberOfPointsInCell(0);
        if (numberOfPointsPerCell == 4)
        {
            invoke(EntropyIndependentGaussian<4, 16>{this->IsoValue}, vtkmDataSet.GetCellSet(), noEnsmeanArray, noEnsstdevArray, noEnscrossProbability, noEnsnumNonZeroProb, noEnsentropy);
        }
        else if (numberOfPointsPerCell == 3)
        {
            invoke(EntropyIndependentGaussian<3, 8>{this->IsoValue}, vtkmDataSet.GetCellSet(), noEnsmeanArray, noEnsstdevArray, noEnscrossProbability, noEnsnumNonZeroProb, noEnsentropy);
        }
        else
        {
            throw std::runtime_error("numberOfPointsPerCell is " + std::to_string(numberOfPointsPerCell) + " is unsupported");
        }
        // only add the entropy array into result
        outputDataSet.AddCellField("entropy_no_ens_" + std::to_string(noEnsId), noEnsentropy);

        // compute the sum of difference compared with the original entropy value
        vtkm::cont::ArrayHandle<vtkm::Float64> diffArray;
        invoke(ComputeDiffSum{}, noEnsentropy, entropy, diffArray);

        outputDataSet.AddCellField("entropy_diff_no_ens_" + std::to_string(noEnsId), diffArray);

        vtkm::Float64 totalDiffValue =
            vtkm::cont::Algorithm::Reduce(diffArray, 0, vtkm::Sum());
        std::cout << totalDiffValue << std::endl;
    }

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
        std::cout << "debug before vtk to vtkm converting" << std::endl;
        // Only the data with the fixed cellset can be converted into the vtkm data set
        vtkm::cont::DataSet in = tovtkm::Convert(input, tovtkm::FieldsFlag::PointsAndCells);
        std::cout << "debug after vtk to vtkm converting" << std::endl;

        std::cout << "debug input vtkm data" << std::endl;
        in.PrintSummary(std::cout);
        vtkm::Id totalNumEnsemble = 0;
        vtkm::cont::DataSet result = this->CallUncertainContourWorklet(in, totalNumEnsemble);

        // Convert the result back from VTKm to VTK.
        // It would be easier if there was a simple method to just convert from general
        // vtkm::cont::DataSet to vtkDataSet. However, that does not exist so you have
        // to copy the vtkDataSet structure in VTK and copy the new fields over. I think
        // it is done this way to prevent creating new arrays for what should be shallow
        // copies. (Maybe in the future the data sharing will be good enough where that
        // is not an issue.)

        // TODO, is there a way to only copy the cellset of the origianl data set?
        // or removing the original ens field in data set
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
        // Copy the field with no_ens value
        for (int noEnsId = 0; noEnsId < totalNumEnsemble; noEnsId++)
        {
            std::string NoEnsName = "entropy_no_ens_" + std::to_string(noEnsId);
            std::string NoEnsDiffName = "entropy_diff_no_ens_" + std::to_string(noEnsId);
            copyField(NoEnsName);
            copyField(NoEnsDiffName);
        }
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