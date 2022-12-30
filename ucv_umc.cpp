#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayCopy.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/DispatcherMapTopology.h>

class ProbMCWorklet : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    ProbMCWorklet(int isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldMinType, typename InPointFieldMaxType, typename OutCellFieldType>

    VTKM_EXEC void operator()(
        const InPointFieldMinType &inPointFieldVecMin,
        const InPointFieldMaxType &inPointFieldVecMax,
        OutCellFieldType &outCellFieldCProb) const
    {
        // how to process the case where there are multiple variables
        // for current data set, there are 8 points for each cell
        // try to make it more general
        vtkm::IdComponent numPoints = inPointFieldVecMin.GetNumberOfComponents();

        vtkm::FloatDefault allPositiveProb = 1.0;
        vtkm::FloatDefault allNegativeProb = 1.0;
        vtkm::FloatDefault allCrossProb = 0.0;

        vtkm::FloatDefault positiveProb;
        vtkm::FloatDefault negativeProb;

        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            vtkm::FloatDefault minV = inPointFieldVecMin[pointIndex];
            vtkm::FloatDefault maxV = inPointFieldVecMax[pointIndex];

            if (this->m_isovalue <= minV)
            {
                positiveProb = 1.0;
                negativeProb = 0.0;
            }
            else if (this->m_isovalue >= maxV)
            {
                positiveProb = 0.0;
                negativeProb = 1.0;
            }
            else
            {
                // assuming we use the uniform distribution
                positiveProb = (maxV - this->m_isovalue) / (maxV - minV);
                negativeProb = 1.0 - positiveProb;
            }

            allPositiveProb *= positiveProb;
            allNegativeProb *= negativeProb;
        }

        allCrossProb = 1 - allPositiveProb - allNegativeProb;

        outCellFieldCProb = allCrossProb;
    }

private:
    int m_isovalue;
};

class ProbMCEntropyWorklet : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    ProbMCEntropyWorklet(int isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldMinType, typename InPointFieldMaxType, typename OutCellFieldType>

    VTKM_EXEC void operator()(
        const InPointFieldMinType &inPointFieldVecMin,
        const InPointFieldMaxType &inPointFieldVecMax,
        OutCellFieldType &outCellFieldCProb) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecMin.GetNumberOfComponents();
        // there are 8 points for each cell

        vtkm::FloatDefault allPositiveProb = 1.0;
        vtkm::FloatDefault allNegativeProb = 1.0;
        vtkm::FloatDefault allCrossProb = 0.0;

        vtkm::FloatDefault positiveProb;
        vtkm::FloatDefault negativeProb;

        std::vector<vtkm::FloatDefault> positiveProbList;
        std::vector<vtkm::FloatDefault> negativeProbList;

        positiveProbList.resize(numPoints);
        negativeProbList.resize(numPoints);

        // TODO, try to make this general
        int totalNumCases = 256;
        std::vector<vtkm::FloatDefault> probHistogram;

        probHistogram.resize(totalNumCases);

        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            vtkm::FloatDefault minV = inPointFieldVecMin[pointIndex];
            vtkm::FloatDefault maxV = inPointFieldVecMax[pointIndex];

            if (this->m_isovalue <= minV)
            {
                positiveProb = 1.0;
                negativeProb = 0.0;
            }
            else if (this->m_isovalue >= maxV)
            {
                positiveProb = 0.0;
                negativeProb = 1.0;
            }
            else
            {
                // assuming we use the uniform distribution
                positiveProb = (maxV - this->m_isovalue) / (maxV - minV);
                negativeProb = 1.0 - positiveProb;
            }

            positiveProbList[pointIndex] = positiveProb;
            negativeProbList[pointIndex] = negativeProb;

            allPositiveProb *= positiveProb;
            allNegativeProb *= negativeProb;
        }

        allCrossProb = 1 - allPositiveProb - allNegativeProb;
        outCellFieldCProb = allCrossProb;
        
        //TODO, how to call the fuction from this worklet
        //traverse(1.0, 0, 0, positiveProbList, negativeProbList, probHistogram);

        //TODO extracting the entropy or other values based on probHistogram
    }

    // using recursive call to go through all possibilities
    void traverse(vtkm::FloatDefault currentProb, int depth, int id,
                  std::vector<vtkm::FloatDefault> &positiveProbList,
                  std::vector<vtkm::FloatDefault> &negativeProbList,
                  std::vector<vtkm::FloatDefault> &probHistogram)
    {
        // TODO, make this as a private variable
        // how to set it as a private variable of the worklet
        if (depth == 8)
        {
            if (id > 256)
            {
                throw std::runtime_error("id is supposed to be 0 to 255");
            }
            probHistogram[id] = currentProb;
            return;
        }
        // two branch for current node
        vtkm::FloatDefault nextPosProb = currentProb * positiveProbList[depth];
        vtkm::FloatDefault nextNegProb = currentProb * negativeProbList[depth];

        traverse(nextPosProb, depth + 1, 1 + (id << 1), positiveProbList, negativeProbList, probHistogram);
        traverse(nextNegProb, depth + 1, id << 1, positiveProbList, negativeProbList, probHistogram);
        return;
    }

private:
    int m_isovalue;
};

using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32>;

// executing the uncertainty marching cube
// based on extracted ensemble data
// refer to this exmaple
// https://gitlab.kitware.com/vtk/vtk-m/-/blob/master/tutorial/point_to_cell.cxx
// using the WorkletVisitCellsWithPoints
int main(int argc, char *argv[])
{
    // init the vtkm (set the backend and log level here)
    vtkm::cont::Initialize(argc, argv);

    if (argc != 4)
    {
        std::cout << "executable <filename> <fieldname> <isovalue>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    int isovalue = std::stoi(argv[3]);
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    std::cout << "check dataset:" << std::endl;
    inData.PrintSummary(std::cout);

    const auto &inFieldRaw = inData.GetField(fieldName);

    vtkm::cont::UnknownArrayHandle outArray;

    using WorkletType = ProbMCWorklet;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;

    DispatcherType dispatcher(ProbMCWorklet{isovalue});

    auto resolveType = [&](const auto &inputArray)
    {
        // try to derive the type of the input field
        using T = typename std::decay_t<decltype(inputArray)>::ValueType;

        // and use the type of raw field as a reminder for predicting derived parameters
        // this shallowIfPossible can cast the variable to dedicated type
        // otherwise, it will use the deep copy to transfer data to dedicated type
        // we use the shallowcopy to fix the type of the input array
        // otherwise, the autonomic type deduction will take lots of time
        // for multiple input field in the worklet

        vtkm::cont::ArrayHandle<T> inFieldMin;
        vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField("ensemble_min").GetData(), inFieldMin);
        vtkm::cont::ArrayHandle<T> inFieldMax;
        vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField("ensemble_max").GetData(), inFieldMax);

        // using the same type as the assumption for the output type
        vtkm::cont::ArrayHandle<T> result;

        dispatcher.Invoke(inData.GetCellSet(), inFieldMin, inFieldMax, result);

        outArray = result;
    };

    inFieldRaw.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);
    // there are some compiling issues for using the CastAndCall direactly
    //vtkm::cont::CastAndCall(inFieldRaw, resolveType);

    std::cout << "===data summary after adding the field array:" << std::endl;
    inData.AddCellField("cross_prob", outArray);
    inData.PrintSummary(std::cout);

    // TODO output the dataset
    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(inData);
    return 0;
}
