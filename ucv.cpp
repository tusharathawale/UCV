#include <vtkm/cont/Initialize.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/cont/ArrayHandle.h>

struct GoThroughNeigoborhood : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    GoThroughNeigoborhood(int neighborhoodSize)
        : m_neighborhoodSize(neighborhoodSize)
    {
    }

    using ControlSignature = void(CellSetIn, FieldInNeighborhood, FieldOut);

    using ExecutionSignature = void(_2,_3);

    // where are these template parameters are specified?
    // is it following the fixed format?
    // the first is the type of input, the second is the type of output

    template <typename InNeighborhoodT, typename T>
    VTKM_EXEC void operator()(const InNeighborhoodT &input,
                              T &output) const
    {
        //how to get the bound of the bounding box of the input data?
        //the type of input is the array ArrayPortalBasicRead
        //we need to skip some of the data points
        
        // the type of input is in src/vtk-m/vtkm/exec/FieldNeighborhood.h
        // std::cout << "this->m_neighborhoodSize " << this->m_neighborhoodSize << std::endl;
        auto flatIndex = input.Get({1,1,1});
        if(flatIndex!=0){
            std::cout << "check the flatindex " << flatIndex << std::endl;
        }
        

        for (vtkm::IdComponent k = -this->m_neighborhoodSize; k <= this->m_neighborhoodSize; ++k)
        {
            for (vtkm::IdComponent j = -this->m_neighborhoodSize; j <=this->m_neighborhoodSize; ++j)
            {
                for (vtkm::IdComponent i = -this->m_neighborhoodSize; i <= this->m_neighborhoodSize; ++i)
                {
                    //is this input change everytime?
                    vtkm::FloatDefault value = static_cast<T>(input.Get(i, j, k));
                    //check the first value of input
                    
                    if(flatIndex!=0 && value!=0){
                        std::cout << "flatIndex " << flatIndex << " test value " << value << std::endl; 
                    }

                    
                }
            }
        }

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
    // the raw_data for testing is 832*832*494
    inData.PrintSummary(std::cout);

    // TODO try to use the point neighborhood
    // skip particular on if it is covered by previous points
    // refer to this:
    // https://gitlab.kitware.com/vtk/vtk-m/-/blob/release-1.9/vtkm/filter/image_processing/worklet/ComputeMoments.h
    // https://gitlab.kitware.com/vtk/vtk-m/-/blob/release-1.9/vtkm/filter/image_processing/ImageMedian.cxx
    // not sure how to use that boundy thing
    /* try the neighborhood worklet */

    std::cout << "get cell set" << std::endl;
    const vtkm::cont::UnknownCellSet &inputCellSet = inData.GetCellSet();
    std::cout << "get field" << std::endl;

    auto field = inData.GetField(fieldName);

    vtkm::cont::UnknownArrayHandle outArray;

    using WorkletType = GoThroughNeigoborhood;
    using DispatcherType = vtkm::worklet::DispatcherPointNeighborhood<WorkletType>;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result;
    std::cout << "dispatcher and invoke" << std::endl;


    auto resolveType = [&](const auto &concrete)
    {
        DispatcherType dispatcher(WorkletType{1});
        dispatcher.Invoke(inputCellSet, concrete, result);
        outArray = result;
    };

    field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    return 0;
}
