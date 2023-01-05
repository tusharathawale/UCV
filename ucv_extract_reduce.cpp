#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletReduceByKey.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/Keys.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <vtkm/worklet/WorkletMapField.h>

#include <float.h>

struct CreateNewKeyWorklet : public vtkm::worklet::WorkletMapField
{

    VTKM_CONT
    CreateNewKeyWorklet(vtkm::Id rawDimx, vtkm::Id rawDimy, vtkm::Id rawDimz,
                        vtkm::Id blocksize) : m_rawDimx(rawDimx), m_rawDimy(rawDimy), m_rawDimz(rawDimz), m_blocksize(blocksize)
    {
        m_numberBlockx = m_rawDimx % blocksize == 0 ? m_rawDimx / blocksize : m_rawDimx / blocksize + 1;
        m_numberBlocky = m_rawDimy % blocksize == 0 ? m_rawDimy / blocksize : m_rawDimy / blocksize + 1;
        m_numberBlockz = m_rawDimz % blocksize == 0 ? m_rawDimz / blocksize : m_rawDimz / blocksize + 1;
    };

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

        // typename OriginalValuesVecType::ComponentType max = 0;

        for (vtkm::IdComponent index = 0;
             index < originalValues.GetNumberOfComponents(); index++)
        {
            vtkm::FloatDefault originalvalue = originalValues[index];
            boxMin = vtkm::Min(boxMin, originalvalue);
            boxMax = vtkm::Max(boxMax, originalvalue);
        }
        // TODO, how to return the min and max here or multiple values
        // TODO, maybe integrating the histogram operations here

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

    if (argc != 3)
    {
        std::cout << "executable <filename> <fieldname>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
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

    // vtkm::cont::ArrayHandle<vtkm::Id> keyArray;
    // keyArray.Allocate(xdim * ydim * zdim);
    // for (vtkm::Id i = 0; i < keyArray.GetNumberOfValues(); ++i)
    //{
    //     keyArray.WritePortal().Set(i, i);
    // }

    auto keyArray =
        vtkm::cont::ArrayHandleCounting<vtkm::Id>(0, 1, static_cast<vtkm::Id>(xdim * ydim * zdim));
    //  the value can be set as a parameter
    vtkm::Id blocksize = 4;
    vtkm::Id numberBlockx = xdim % blocksize == 0 ? xdim / blocksize : xdim / blocksize + 1;
    vtkm::Id numberBlocky = ydim % blocksize == 0 ? ydim / blocksize : ydim / blocksize + 1;
    vtkm::Id numberBlockz = zdim % blocksize == 0 ? zdim / blocksize : zdim / blocksize + 1;
    /* the sequential way for creating the key array
        std::cout << "number blocks " << numberBlockx << " " << numberBlocky << " " << numberBlockz << std::endl;

        // TODO computing number of coarse grid at each dim
        // computing the id of the coarse grid in the subsequent
        // for loops

        vtkm::Id index = 0;
        for (vtkm::Id z = 0; z < zdim; z++)
        {
            vtkm::Id newidz = z / blocksize;
            for (vtkm::Id y = 0; y < ydim; y++)
            {
                vtkm::Id newidy = y / blocksize;
                for (vtkm::Id x = 0; x < xdim; x++)
                {
                    vtkm::Id newidx = x / blocksize;

                    vtkm::Id key = newidx + newidy * numberBlockx + newidz * numberBlockx * numberBlocky;
                    //  set the key to array
                    keyArray.WritePortal().Set(index, key);
                    index++;
                }
            }
        }
    */
    vtkm::cont::ArrayHandle<vtkm::Id> keyArrayNew;
    using DispatcherCreateKey = vtkm::worklet::DispatcherMapField<CreateNewKeyWorklet>;
    DispatcherCreateKey dispatcher(CreateNewKeyWorklet{xdim, ydim, zdim, blocksize});

    // auto resolveType = [&](const auto &concrete)
    //{
    dispatcher.Invoke(keyArray, keyArrayNew);
    //};

    // keyArray.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
    //     resolveType);

    // Try to use another worklet to generate this key array

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
    reducedDataSet.AddPointField("ensemble_min", minArray);
    reducedDataSet.AddPointField("ensemble_max", maxArray);
    reducedDataSet.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + std::string("_ReduceDerived.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(reducedDataSet);

    return 0;
}
