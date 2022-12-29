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
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell);

    // using ControlSignature = void(CellSetIn,
    //                               FieldInPoint,
    //                               FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4, _5);
    // using ExecutionSignature = void(_2, _3);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldType1, typename InPointFieldType2, typename InPointFieldType3, typename OutCellFieldType>
    // VTKM_EXEC void operator()(const InPointFieldType &inPointFieldVecRaw,
    //                           OutCellFieldType &outCellField) const
    VTKM_EXEC void operator()(const InPointFieldType1 &inPointFieldVecRaw,
                              const InPointFieldType2 &inPointFieldVecMin,
                              const InPointFieldType3 &inPointFieldVecMax,
                              OutCellFieldType &outCellField) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecRaw.GetNumberOfComponents();

        // std::cout << "numPoints for each input" << numPoints << std::endl;

        outCellField = OutCellFieldType(0);

        // how to go through different field?
        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            // outCellField = outCellField + inPointFieldVec[pointIndex];
            // TODO operations about marching cube things
            // vtkm::FloatDefault raw = inPointFieldVec[pointIndex][0];
            // vtkm::FloatDefault fieldMin = inPointFieldVec[pointIndex][1];
            // vtkm::FloatDefault fieldMax = inPointFieldVec[pointIndex][2];
        }

        outCellField =
            static_cast<OutCellFieldType>(outCellField / static_cast<vtkm::FloatDefault>(numPoints));
    }

private:
    int m_isovalue;
};

// executing the uncertainty marching cube
// based on extracted ensemble data
// refer to this exmaple
// https://gitlab.kitware.com/vtk/vtk-m/-/blob/master/tutorial/point_to_cell.cxx
// using the WorkletVisitCellsWithPoints
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
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    std::cout << "check dataset:" << std::endl;
    inData.PrintSummary(std::cout);

    const auto &inFieldRaw = inData.GetField(fieldName);
    const auto &inFieldMin = inData.GetField("min");
    const auto &inFieldMax = inData.GetField("max");

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> outArray;

    using WorkletType = ProbMCWorklet;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;

    DispatcherType dispatcher(ProbMCWorklet{900});

    // how to process the case where there are multiple fields
    // maybe gradient filter is an example
    // we currently do not use the resolveType, just call the Invoke direactly with concrete type?

    // dispatcher.Invoke(inData.GetCellSet(), inFieldRaw, inFieldMin, inFieldMax, outArray);
    // it takes comparatively long time to compile this
    // around several minutes for doing this, not sure the reason
    // dispatcher.Invoke(inData.GetCellSet(), inFieldRaw, outArray);
    dispatcher.Invoke(inData.GetCellSet(), inFieldRaw, inFieldMin, inFieldMax, outArray);

    // TODO add results into the dataset

    // TODO output the dataset

    return 0;
}
