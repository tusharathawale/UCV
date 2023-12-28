#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletMapField.h>
#include "ucvworklet/EntropyUniform.hpp"

struct AmrexMinMax : public vtkm::worklet::WorkletMapField
{
    AmrexMinMax(const vtkm::FloatDefault error) : Error(error){};

    using ControlSignature = void(FieldIn, FieldOut, FieldOut);
    using ExecutionSignature = void(_1, _2, _3);

    template <typename EnputType, typename OutputType>
    VTKM_EXEC void operator()(
        const EnputType &originalValue, OutputType &minValue, OutputType &maxValue) const
    {
        minValue = originalValue - this->Error;
        maxValue = originalValue + this->Error;
    }

    vtkm::FloatDefault Error = 0;
};

int main(int argc, char *argv[])
{

    // init the vtkm (set the backend and log level here)
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    vtkm::cont::Timer timer{initResult.Device};
    std::cout << "initResult.Device: " << initResult.Device.GetName() << " timer device: " << timer.GetDevice().GetName() << std::endl;

    if (argc != 5)
    {
        std::cout << "executable [VTK-m options] <filename> <fieldname> <isovalue> <error value>" << std::endl;
        std::cout << "VTK-m options are:\n";
        std::cout << initResult.Usage << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    double isovalue = std::atof(argv[3]);
    double errvalue = std::atof(argv[4]);

    std::cout << "check input fileName " << fileName << "fieldName " << fieldName << " isovalue " << isovalue << " errvalue " << errvalue << std::endl;

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();
    inData.PrintSummary(std::cout);

    vtkm::Id dx = 128;
    vtkm::Id dy = 128;
    vtkm::Id dz = 1024;

    // get min and max through error data
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> minArray;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxArray;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> entropyResult;

    vtkm::cont::Invoker invoke;
    auto resolveType = [&](auto &concreteArray)
    {
        invoke(AmrexMinMax{isovalue}, concreteArray, minArray, maxArray);
    };

    inData.GetField(fieldName)
        .GetData()
        .CastAndCallWithExtractedArray(resolveType);

    // call the min and max filter
    invoke(EntropyUniform{isovalue}, inData.GetCellSet(), minArray, maxArray, crossProb, numNonZeroProb, entropyResult);

    // output data
    vtkm::cont::DataSetBuilderUniform dataSetBuilder;
    const vtkm::Id3 reducedDims(dx, dy, dz);
    const vtkm::Id3 reducedOrigin(0, 0, 0);
    const vtkm::Id3 reducedSpaceing(1, 1, 1);
    vtkm::cont::DataSet outputDataSet = dataSetBuilder.Create(reducedDims, reducedOrigin, reducedSpaceing);

    outputDataSet.AddCellField("cross_prob", crossProb);
    outputDataSet.AddCellField("num_nonzero_prob", numNonZeroProb);
    outputDataSet.AddCellField("entropy", entropyResult);
    outputDataSet.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    std::string outputFileName = "amrex_output.vtk";
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(outputDataSet);

    return 0;
}