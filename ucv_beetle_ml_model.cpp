#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletMapField.h>

#include "ucvworklet/EntropyUniform.hpp"
#include "ucvworklet/EntropyIndependentGaussian.hpp"

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayCopy.h>

struct HurricanMinMax : public vtkm::worklet::WorkletMapField
{
    HurricanMinMax(const vtkm::FloatDefault error) : Error(error){};

    using ControlSignature = void(FieldIn, FieldOut, FieldOut);
    using ExecutionSignature = void(_1, _2, _3);

    template <typename EnputType, typename OutputType>
    VTKM_EXEC void operator()(
        const EnputType &originalValue, OutputType &minValue, OutputType &maxValue) const
    {
        minValue = originalValue - this->Error / 2.0;
        maxValue = originalValue + this->Error / 2.0;
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

    std::string fileName = "/Users/zhewang/Downloads/recon_beetle.vtk";
    std::string fieldNameAvg = "Average";
    std::string fieldNameStdv = "Standard%20Deviation";
    vtkm::FloatDefault isovalue = 1000;

    std::cout << "check input fileName " << fileName << "fieldNameAvg " << fieldNameAvg << " fieldNameStdv " << fieldNameStdv << std::endl;

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();
    inData.PrintSummary(std::cout);

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> entropyResult;

    vtkm::cont::ArrayHandle<vtkm::Float64> stdDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(fieldNameStdv).GetData(), stdDataArray);

    // using indepedent gaussian distribution
    auto resolveTypeIG = [&](const auto &concrete)
    {
        using WorkletType = EntropyIndependentGaussian<8, 256>;
        using DispatcherEntropyIG = vtkm::worklet::DispatcherMapTopology<WorkletType>;

        DispatcherEntropyIG dispatcherEntropyIG(EntropyIndependentGaussian<8, 256>{isovalue});
        // assume compressed data is mean
        dispatcherEntropyIG.Invoke(inData.GetCellSet(), concrete, stdDataArray, crossProb, numNonZeroProb, entropyResult);
    };

    inData.GetField(fieldNameAvg).GetData().CastAndCallWithExtractedArray(resolveTypeIG);

    // output data
    vtkm::Id dx = 247;
    vtkm::Id dy = 416;
    vtkm::Id dz = 416;
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
    std::string outputFileName = "beetle_ml_output_uni_iso_" + std::to_string(isovalue) + ".vtk";
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(outputDataSet);

    return 0;
}