#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandle.h>

#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/DispatcherMapTopology.h>

class ProbMCWorklet : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    ProbMCWorklet(int isovalue)
        : m_isovalue(isovalue){};

    // two input variables
    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4, _5);
    // the first parameter is binded with the worklet
    using InputDomain = _1;

    // Gradient filter is an example that recieves multiple fields
    template <typename InPointFieldType, typename InPointFieldMinType, typename InPointFieldMaxType, typename OutCellFieldType>
    VTKM_EXEC void operator()(const InPointFieldType &inPointFieldRaw,
                              const InPointFieldMinType &inPointFieldMin,
                              const InPointFieldMaxType &inPointFieldMax,
                              OutCellFieldType &outCellField) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldRaw.GetNumberOfComponents();

        std::cout << "numPoints for each input" << numPoints << std::endl;

        outCellField = OutCellFieldType(0);

        // howt to go through each point in the cell?
        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            // outCellField = outCellField + inPointFieldVec[pointIndex];
            // TODO marching cube things
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
    auto &inFieldRaw = inData.GetField(fieldName);
    auto &inFieldMin = inData.GetField("min");
    auto &inFieldMax = inData.GetField("max");

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
    dispatcher.Invoke(inData.GetCellSet(), inFieldRaw, inFieldMin, inFieldMax, outArray);

    // TODO add results into the dataset

    // TODO output the dataset

    return 0;
}
