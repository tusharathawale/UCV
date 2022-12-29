// this is the depricated code, it uses vector to transfer data which is not good 
// use this as a reference
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandle.h>

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
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3);
    // using ExecutionSignature = void(_2, _3);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType2 is supposed to be a vector of vector
    // the first level vector is 8 points
    // the second level vector represents derived values
    template <typename InPointFieldType, typename OutCellFieldType>
    VTKM_EXEC void operator()(const InPointFieldType &inPointFieldVecDerived,
                              OutCellFieldType &outCellField) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecDerived.GetNumberOfComponents();
        // there are 8 points for each cell

        vtkm::FloatDefault allPositiveProb = 1.0;
        vtkm::FloatDefault allNegativeProb = 1.0;
        vtkm::FloatDefault allCrossProb = 0.0;

        vtkm::FloatDefault positiveProb;
        vtkm::FloatDefault negativeProb;

        vtkm::FloatDefault entropy;
        vtkm::FloatDefault entropyLocal;

        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            vtkm::Vec<vtkm::FloatDefault, 2> derivedVec = inPointFieldVecDerived[pointIndex];
            vtkm::FloatDefault minV = derivedVec[0];
            vtkm::FloatDefault maxV = derivedVec[1];

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

            // for each vertex, there are two probabilities
            // [pprob, nprob]
        }

        //create 256 cases
        //it requires a process to go through a binary tree with the depth equals to 9
        //based on the information on each leave node 
        //extract entropy or other parameters, such as number of possible cases
        //there should be a worklet that create this histogram
        //for each cell, there are 256 possible values
        //some possible trim operation:
        //trim the branch with the prob=0

        allCrossProb = 1 - allPositiveProb - allNegativeProb;

        outCellField[0] = allPositiveProb;
        outCellField[1] = allNegativeProb;
        outCellField[2] = allCrossProb;

        // how to compute histogram based on the case table
        // computing the histogram? the key of case histogram is the case, value is its probability
        // There are 256 cases per cell
        // refering to the marching cube worklet
        // and value such as entropy
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
                                  vtkm::UInt32,
                                  vtkm::Vec<vtkm::FloatDefault, 2>>;

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
    
    // it is potential to trigger long compiling time if the type is not clear
    const auto &inFieldRaw = inData.GetField(fieldName);
    // vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 2>> inFieldDerived = inData.GetField("derived").GetData();
    const auto &inFieldDerived = inData.GetField("derived");

    std::cout << "ok to get field, call dispatcher" << std::endl;

    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 4>> outArray;

    using WorkletType = ProbMCWorklet;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;

    DispatcherType dispatcher(ProbMCWorklet{isovalue});

    // how to process the case where there are multiple fields
    // maybe gradient filter is an example
    // we currently do not use the resolveType, just call the Invoke direactly with concrete type?

    // dispatcher.Invoke(inData.GetCellSet(), inFieldRaw, inFieldMin, inFieldMax, outArray);
    // it takes comparatively long time to compile this
    // around several minutes for doing this, not sure the reason
    // dispatcher.Invoke(inData.GetCellSet(), inFieldRaw, outArray);

    // question, how to use the field callback operation if there are multiple field

    auto resolveType = [&](const auto &concreteDerived)
    {
        dispatcher.Invoke(inData.GetCellSet(), concreteDerived, outArray);
    };

    // we do not need the raw value to compute the histogram things

    inFieldDerived.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    std::cout << "data summary after adding the field array:" << std::endl;
    inData.AddCellField("prob", outArray);
    inData.PrintSummary(std::cout);

    // TODO output the dataset
    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(inData);
    return 0;
}
