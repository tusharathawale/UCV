#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/cont/ArrayHandle.h>

#include <vtkm/cont/DataSetBuilderUniform.h>

#include <float.h>

struct GoThroughNeigoborhood : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    GoThroughNeigoborhood(int neighborhoodSize)
        : m_neighborhoodSize(neighborhoodSize){};

    using ControlSignature = void(CellSetIn, FieldInNeighborhood, FieldOut, FieldOut);

    using ExecutionSignature = void(_2, _3, _4);

    template <typename InNeighborhoodT, typename T>
    VTKM_EXEC void operator()(const InNeighborhoodT &input,
                              T &output1, T &output2) const
    {

        // some questions:
        // 1. how to process the case where the ijk is out of the box
        // 2. how to get the coarse results? Large data -> small data
        // 3. how to combine it with the ascent that contains multiple blocks?
        vtkm::FloatDefault boxMin = DBL_MAX;
        // the min value for the field is 0 in this dataset
        // we can reset this number if it is different
        vtkm::FloatDefault boxMax = 0;

        for (vtkm::IdComponent k = -this->m_neighborhoodSize; k <= this->m_neighborhoodSize; ++k)
        {
            for (vtkm::IdComponent j = -this->m_neighborhoodSize; j <= this->m_neighborhoodSize; ++j)
            {
                for (vtkm::IdComponent i = -this->m_neighborhoodSize; i <= this->m_neighborhoodSize; ++i)
                {
                    // the type of T1 and T2 are same
                    vtkm::FloatDefault value = static_cast<T>(input.Get(i, j, k));

                    // add other values as needed for the next step
                    boxMax = std::max(boxMax, value);
                    boxMin = std::min(boxMin, value);
                }
            }
        }
        output1 = boxMin;
        output2 = boxMax;
    }

private:
    int m_neighborhoodSize;
};

using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32>;

// TODO, extracting the ensemble data based on vtk operation
int main(int argc, char *argv[])
{
    // init the vtkm (set the backend and log level here)
    vtkm::cont::Initialize(argc, argv);

    if (argc != 4)
    {
        std::cout << "executable <filename> <fieldname1> <fieldname2>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName1 = argv[2];
    std::string fieldName2 = argv[3];
    // load the dataset (beetles data set, structured one)
    // TODO, the data set can be distributed between different ranks

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // check the property of the data
    inData.PrintSummary(std::cout);

    // std::cout << "get cell set" << std::endl;
    const vtkm::cont::UnknownCellSet &inputCellSet = inData.GetCellSet();
    // std::cout << "get field" << std::endl;

    auto field1 = inData.GetField(fieldName1);
    auto field2 = inData.GetField(fieldName2);

    using WorkletType = GoThroughNeigoborhood;
    using DispatcherType = vtkm::worklet::DispatcherPointNeighborhood<WorkletType>;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result11;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result12;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result21;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result22;

    auto resolveType1 = [&](const auto &concrete)
    {
        DispatcherType dispatcher(WorkletType{2});
        dispatcher.Invoke(inputCellSet, concrete, result11, result12);
    };

    field1.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType1);

    auto resolveType2 = [&](const auto &concrete)
    {
        DispatcherType dispatcher(WorkletType{2});
        dispatcher.Invoke(inputCellSet, concrete, result21, result22);
    };

    field2.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType2);

    // add new array into data set
    std::cout << "data summary after adding the field array:" << std::endl;

    // it is not good to use vector to store multiple data
    // the code can be bespoke
    // we can add more derived value here
    // such as mean and stdev used for gaussian distribution

    // creating new data set
    auto cellSet = inData.GetCellSet();

    vtkm::cont::CellSetStructured<3> structCellSet =
        cellSet.AsCellSet<vtkm::cont::CellSetStructured<3>>();

    vtkm::Id3 pointDims = structCellSet.GetPointDimensions();

    std::cout << "------" << std::endl;
    std::cout << "point dim: " << pointDims[0] << " " << pointDims[1] << " " << pointDims[2] << std::endl;

    // go through all points and set the specific key
    vtkm::Id xdim = pointDims[0];
    vtkm::Id ydim = pointDims[1];
    vtkm::Id zdim = pointDims[2];
    const vtkm::Id3 dims(xdim, ydim, zdim);
    vtkm::Bounds bounds(-1.02655e+09, 9.42128e+08, -1.02353e+09, 1.00073e+09, -9.16584e+08, 9.32988e+08);
    auto reducedOrigin = bounds.MinCorner();
    vtkm::FloatDefault spacex = (bounds.X.Max - bounds.X.Min) / (xdim - 1);
    vtkm::FloatDefault spacey = (bounds.Y.Max - bounds.Y.Min) / (ydim - 1);
    vtkm::FloatDefault spacez = (bounds.Z.Max - bounds.Z.Min) / (zdim - 1);

    vtkm::Vec3f_64 reducedSpaceing(spacex, spacey, spacez);

    vtkm::cont::DataSetBuilderUniform dataSetBuilder;

    vtkm::cont::DataSet datasetOutput = dataSetBuilder.Create(dims, reducedOrigin, reducedSpaceing);

    datasetOutput.AddPointField(fieldName1 + "_ensemble_min", result11);
    datasetOutput.AddPointField(fieldName1 + "_ensemble_max", result12);

    datasetOutput.AddPointField(fieldName2 + "_ensemble_min", result21);
    datasetOutput.AddPointField(fieldName2 + "_ensemble_max", result22);

    datasetOutput.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + std::string("_ReduceMinMax_2Variables.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(datasetOutput);

    return 0;
}
