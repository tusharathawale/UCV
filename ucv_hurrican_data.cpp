#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletMapField.h>

#include "ucvworklet/EntropyUniform.hpp"
#include "ucvworklet/EntropyIndependentGaussian.hpp"

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

    if (argc != 6)
    {
        std::cout << "executable [VTK-m options] <filename> <fieldname> <isovalue> <error value> <std value>" << std::endl;
        std::cout << "VTK-m options are:\n";
        std::cout << initResult.Usage << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    double isovalue = std::atof(argv[3]);
    double errvalue = std::atof(argv[4]);
    double stdvalue = std::atof(argv[5]);

    std::cout << "check input fileName " << fileName << "fieldName " << fieldName << " isovalue " << isovalue << " errvalue " << errvalue  << " stdev " << stdvalue << std::endl;

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();
    inData.PrintSummary(std::cout);

    vtkm::Id dx = 496;
    vtkm::Id dy = 496;
    vtkm::Id dz = 96;

    // get min and max through error data
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> minArray;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxArray;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> entropyResult;

    vtkm::cont::Invoker invoke;
    auto resolveType = [&](auto &concreteArray)
    {
        invoke(HurricanMinMax{errvalue}, concreteArray, minArray, maxArray);
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
    std::string outputFileName = "hurrican_output_uni_iso" + std::to_string(isovalue) + "_err" + std::to_string((errvalue)) + ".vtk";
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(outputDataSet);

    // reset the buffer
    crossProb.Fill(0);
    numNonZeroProb.Fill(0);
    entropyResult.Fill(0);

    // using indepedent gaussian distribution, mean value is same with the compressed value
    std::cout << "------computing gaussian output" << std::endl;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdevArray;
    double stdev = stdvalue;

    stdevArray.AllocateAndFill(dx * dy * dz, stdev);

    auto resolveTypeIG = [&](const auto &concrete)
    {
        using WorkletType = EntropyIndependentGaussian<8, 256>;
        using DispatcherEntropyIG = vtkm::worklet::DispatcherMapTopology<WorkletType>;

        DispatcherEntropyIG dispatcherEntropyIG(EntropyIndependentGaussian<8, 256>{isovalue});
        // assume compressed data is mean
        dispatcherEntropyIG.Invoke(inData.GetCellSet(), concrete, stdevArray, crossProb, numNonZeroProb, entropyResult);
    };

    inData.GetField(fieldName).GetData().CastAndCallWithExtractedArray(resolveTypeIG);

    vtkm::cont::DataSet outputDataSetGaussian = dataSetBuilder.Create(reducedDims, reducedOrigin, reducedSpaceing);

    outputDataSetGaussian.AddCellField("cross_prob", crossProb);
    outputDataSetGaussian.AddCellField("num_nonzero_prob", numNonZeroProb);
    outputDataSetGaussian.AddCellField("entropy", entropyResult);
    outputDataSetGaussian.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    std::string outputFileNameGS = "hurrican_output_gaussian_iso" + std::to_string(isovalue) + "_std" + std::to_string((stdev)) + ".vtk";
    vtkm::io::VTKDataSetWriter writegs(outputFileNameGS);
    writegs.WriteDataSet(outputDataSetGaussian);

    return 0;
}

/*
for the hurricane_compressed
The data is in 500*500*100, fp32. And the iso-value I used is 0.0006.
The error bound is 0.00046 and the Variance is 1.15E-10 (the standard deviation is 0.0000107).  0.000057 for old data
For compressed-larger-eb data, the stdev is 0.000097
*/