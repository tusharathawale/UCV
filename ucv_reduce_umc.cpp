#include <float.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletReduceByKey.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/Keys.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/DispatcherMapTopology.h>

class ProbMCEntropyWorklet : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    ProbMCEntropyWorklet(int isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4, _5, _6);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldMinType, typename InPointFieldMaxType, typename OutCellFieldType1, typename OutCellFieldType2, typename OutCellFieldType3>

    VTKM_EXEC void operator()(
        const InPointFieldMinType &inPointFieldVecMin,
        const InPointFieldMaxType &inPointFieldVecMax,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy) const
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

        // there are 2^n total cases
        int totalNumCases = static_cast<int>(vtkm::Pow(2.0, static_cast<vtkm::FloatDefault>(numPoints)));
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

        // TODO, use the number of vertesies as another parameter
        traverse(1.0, 0, 0, numPoints, positiveProbList, negativeProbList, probHistogram);

        // extracting the entropy or other values based on probHistogram

        vtkm::FloatDefault entropyValue = 0;
        vtkm::Id nonzeroCases = 0;
        vtkm::FloatDefault templog = 0;
        for (int i = 0; i < totalNumCases; i++)
        {
            templog = 0;
            if (probHistogram[i] > 0.00001)
            {
                nonzeroCases++;
                templog = vtkm::Log2(probHistogram[i]);
            }
            entropyValue = entropyValue + (-probHistogram[i]) * templog;
        }
        outCellFieldNumNonzeroProb = nonzeroCases;
        outCellFieldEntropy = entropyValue;
    }

    // using recursive call to go through all possibilities
    void traverse(vtkm::FloatDefault currentProb, int depth, int id, const int numPoints,
                  std::vector<vtkm::FloatDefault> &positiveProbList,
                  std::vector<vtkm::FloatDefault> &negativeProbList,
                  std::vector<vtkm::FloatDefault> &probHistogram) const
    {
        // TODO, make this as a private variable
        // how to set it as a private variable of the worklet
        if (depth == numPoints)
        {
            // if (id > 256) how to set this as a worklet parameter?
            //{
            //     throw std::runtime_error("id is supposed to be 0 to 255");
            // }
            probHistogram[id] = currentProb;
            return;
        }
        // two branch for current node
        vtkm::FloatDefault nextPosProb = currentProb * positiveProbList[depth];
        vtkm::FloatDefault nextNegProb = currentProb * negativeProbList[depth];

        traverse(nextPosProb, depth + 1, 1 + (id << 1), numPoints, positiveProbList, negativeProbList, probHistogram);
        traverse(nextNegProb, depth + 1, id << 1, numPoints, positiveProbList, negativeProbList, probHistogram);
        return;
    }

private:
    int m_isovalue;
};

struct CreateNewKeyWorklet : public vtkm::worklet::WorkletMapField
{

    VTKM_CONT
    CreateNewKeyWorklet(vtkm::Id rawDimx, vtkm::Id rawDimy, vtkm::Id rawDimz,
                        vtkm::Id numberBlockx, vtkm::Id numberBlocky, vtkm::Id numberBlockz,
                        vtkm::Id blocksize) : m_rawDimx(rawDimx), m_rawDimy(rawDimy), m_rawDimz(rawDimz),
                                              m_numberBlockx(numberBlockx), m_numberBlocky(numberBlocky), m_numberBlockz(numberBlockz),
                                              m_blocksize(blocksize){};

    vtkm::Id m_rawDimx;
    vtkm::Id m_rawDimy;
    vtkm::Id m_rawDimz;

    vtkm::Id m_numberBlockx;
    vtkm::Id m_numberBlocky;
    vtkm::Id m_numberBlockz;

    vtkm::Id m_blocksize;

    typedef void ControlSignature(FieldIn, FieldOut);
    typedef void ExecutionSignature(_1, _2);

    template <typename T>
    VTKM_EXEC void operator()(const T &globalIndex, vtkm::Id &hixelId) const
    {
        vtkm::Id globalId = static_cast<vtkm::Id>(globalIndex);
        // change globalId to the hixel index
        // compute the raw index firstly
        vtkm::Id rawx = globalId % m_rawDimx;
        vtkm::Id rawy = (globalId / m_rawDimx) % m_rawDimy;
        vtkm::Id rawz = globalId / m_rawDimx / m_rawDimy;

        vtkm::Id newidx = rawx / m_blocksize;
        vtkm::Id newidy = rawy / m_blocksize;
        vtkm::Id newidz = rawz / m_blocksize;

        hixelId = newidx + newidy * m_numberBlockx + newidz * m_numberBlockx * m_numberBlocky;
    };
};

struct ExtractingDistribution : public vtkm::worklet::WorkletReduceByKey
{
    using ControlSignature = void(KeysIn, ValuesIn, ReducedValuesOut, ReducedValuesOut);
    using ExecutionSignature = void(_2, _3, _4);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputType &minValue, OutputType &maxValue) const
    {
        vtkm::FloatDefault boxMin = DBL_MAX;
        vtkm::FloatDefault boxMax = 0;

        for (vtkm::IdComponent index = 0;
             index < originalValues.GetNumberOfComponents(); index++)
        {
            vtkm::FloatDefault originalvalue = originalValues[index];
            boxMin = vtkm::Min(boxMin, originalvalue);
            boxMax = vtkm::Max(boxMax, originalvalue);
        }

        minValue = boxMin;
        maxValue = boxMax;
    }
};

using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32,
                                  vtkm::Id>;

int main(int argc, char *argv[])
{
    // init the vtkm (set the backend and log level here)
    vtkm::cont::Initialize(argc, argv);

    if (argc != 5)
    {
        std::cout << "executable <filename> <fieldname> <blocksize> <isovalue>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    int blocksize = std::stoi(argv[3]);
    int isovalue = std::stoi(argv[4]);

    // load the dataset (beetles data set, structured one)
    // TODO, the data set can be distributed between different ranks

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // check the property of the data
    inData.PrintSummary(std::cout);

    auto field = inData.GetField(fieldName);

    auto cellSet = inData.GetCellSet();

    // Assuming the imput data is the structured data

    bool isStructured = cellSet.IsType<vtkm::cont::CellSetStructured<3>>();
    if (!isStructured)
    {
        std::cout << "the extraction only works for CellSetStructured<3>" << std::endl;
        exit(0);
    }

    vtkm::cont::CellSetStructured<3> structCellSet =
        cellSet.AsCellSet<vtkm::cont::CellSetStructured<3>>();

    vtkm::Id3 pointDims = structCellSet.GetPointDimensions();

    std::cout << "------" << std::endl;
    std::cout << "point dim: " << pointDims[0] << " " << pointDims[1] << " " << pointDims[2] << std::endl;

    // go through all points and set the specific key
    vtkm::Id xdim = pointDims[0];
    vtkm::Id ydim = pointDims[1];
    vtkm::Id zdim = pointDims[2];

    auto keyArray =
        vtkm::cont::ArrayHandleCounting<vtkm::Id>(0, 1, static_cast<vtkm::Id>(xdim * ydim * zdim));

    vtkm::Id numberBlockx = xdim % blocksize == 0 ? xdim / blocksize : xdim / blocksize + 1;
    vtkm::Id numberBlocky = ydim % blocksize == 0 ? ydim / blocksize : ydim / blocksize + 1;
    vtkm::Id numberBlockz = zdim % blocksize == 0 ? zdim / blocksize : zdim / blocksize + 1;

    vtkm::cont::ArrayHandle<vtkm::Id> keyArrayNew;

    // Try to use another worklet to generate this key array
    using DispatcherCreateKey = vtkm::worklet::DispatcherMapField<CreateNewKeyWorklet>;
    DispatcherCreateKey dispatcher(CreateNewKeyWorklet{xdim, ydim, zdim,
                                                       numberBlockx, numberBlocky, numberBlockz,
                                                       blocksize});

    dispatcher.Invoke(keyArray, keyArrayNew);

    // add key array into the dataset, and check the output
    // inData.AddPointField("keyArray", keyArrayNew);
    // std::cout << "------" << std::endl;
    // inData.PrintSummary(std::cout);

    // std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    // std::string outputFileName = fileSuffix + std::string("_Key.vtk");
    // vtkm::io::VTKDataSetWriter write(outputFileName);
    // write.WriteDataSet(inData);

    // try to use the reduce by key to extract the key information of each small blocks
    using DispatcherType = vtkm::worklet::DispatcherReduceByKey<ExtractingDistribution>;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> minArray;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxArray;
    vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

    auto resolveType = [&](const auto &concrete)
    {
        DispatcherType dispatcher;
        dispatcher.Invoke(keys, concrete, minArray, maxArray);
    };

    field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    // create the new data sets for the reduced data
    // the dims for new data sets are numberBlockx*numberBlocky*numberBlockz
    const vtkm::Id3 reducedDims(numberBlockx, numberBlocky, numberBlockz);

    vtkm::cont::DataSetBuilderUniform dataSetBuilder;
    vtkm::cont::DataSet reducedDataSet = dataSetBuilder.Create(reducedDims);

    // generate the new data sets with min and max
    // reducedDataSet.AddPointField("ensemble_min", minArray);
    // reducedDataSet.AddPointField("ensemble_max", maxArray);
    // reducedDataSet.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    // std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    // std::string outputFileName = fileSuffix + std::string("_ReduceDerived.vtk");
    // vtkm::io::VTKDataSetWriter write(outputFileName);
    // write.WriteDataSet(reducedDataSet);

    // we do not write out the dataset, try to call the umc here for the new dataset

    using WorkletType = ProbMCEntropyWorklet;
    using DispatcherProbMCType = vtkm::worklet::DispatcherMapTopology<WorkletType>;

    DispatcherProbMCType dispatcherProbMC(ProbMCEntropyWorklet{isovalue});

    // using the same type as the assumption for the output type
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result1;
    vtkm::cont::ArrayHandle<vtkm::Id> result2;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result3;

    dispatcherProbMC.Invoke(reducedDataSet.GetCellSet(), minArray, maxArray, result1, result2, result3);

    std::cout << "===data summary after adding the field array:" << std::endl;
    reducedDataSet.AddCellField("cross_prob", result1);
    reducedDataSet.AddCellField("num_nonzero_prob", result2);
    reducedDataSet.AddCellField("entropy", result3);

    reducedDataSet.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(reducedDataSet);

    return 0;
}
