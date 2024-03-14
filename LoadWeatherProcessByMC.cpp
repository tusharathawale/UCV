#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>

#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandleRuntimeVec.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleRuntimeVec.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Algorithm.h>

#include "ucvworklet/ExtractMinMaxByErr.hpp"
#include "ucvworklet/CriticalPointWorklet.hpp"
#include "ucvworklet/CriticalPointMonteCarlo.hpp"

using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32,
                                  vtkm::Id>;

void callCriticalPointWorklet(vtkm::cont::DataSet &vtkmDataSet, vtkm::Float64 err, vtkm::Id numSamples)
{

    auto concreteArray = vtkmDataSet.GetField("TestField")
                             .GetData()
                             .AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::Float32>>();

    vtkm::cont::Invoker invoke;
    // Get min and max for each point
    vtkm::cont::ArrayHandle<vtkm::Float64> fieldMin;
    vtkm::cont::ArrayHandle<vtkm::Float64> fieldMax;
    invoke(ExtractMinMaxByErr{err}, concreteArray, fieldMin, fieldMax);

    printSummary_ArrayHandle(fieldMin, std::cout, false);
    printSummary_ArrayHandle(fieldMax, std::cout, false);

    vtkmDataSet.AddPointField("fieldMin", fieldMin);
    vtkmDataSet.AddPointField("fieldMax", fieldMax);

    vtkm::cont::ArrayHandle<vtkm::Float64> outMinProb;
    // Use point neighborhood to go through data
    // invoke(CriticalPointWorklet{}, vtkmDataSet.GetCellSet(), fieldMin, fieldMax, outMinProb);
    invoke(CriticalPointMonteCarloWorklet{numSamples}, vtkmDataSet.GetCellSet(), fieldMin, fieldMax, outMinProb);
    // std::cout << "debug outMinProb:" << std::endl;
    // printSummary_ArrayHandle(outMinProb, std::cout, true);
    vtkmDataSet.AddPointField("MinProb", outMinProb);

    return;
}

int main(int argc, char *argv[])
{
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    vtkm::cont::Timer timer{initResult.Device};

    if (argc != 8)
    {
        //./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
        std::cout << "<executable> <filename> <FieldName> <Dimx> <Dimy> <Dimz> <err> <num samples>" << std::endl;
        exit(0);
    }

    std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

    std::string fileName = std::string(argv[1]);
    std::string fieldName = std::string(argv[2]);

    int dimx = std::stoi(argv[3]);
    int dimy = std::stoi(argv[4]);
    int dimz = std::stoi(argv[5]);

    vtkm::Float64 err = std::stod(argv[6]);
    vtkm::Id numSamples = std::stod(argv[7]);

    const vtkm::Id3 dims(dimx, dimy, dimz);
    // vtkm::cont::DataSetBuilderUniform dataSetBuilder;
    // vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);

    // load input
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // inData.PrintSummary(std::cout);

    // get the global min and max
    auto inTestArray = inData.GetField("TestField").GetData().AsArrayHandle<vtkm::cont::ArrayHandle<vtkm::Float32>>();
    auto gloablMaxValue = vtkm::cont::Algorithm::Reduce(inTestArray, 0, vtkm::Maximum());
    auto gloablMinValue = vtkm::cont::Algorithm::Reduce(inTestArray, vtkm::Infinity64(), vtkm::Minimum());
    // the absolute eb = 1e-2 * (maxv- minv)
    std::cout << "gloablMinValue is " << gloablMinValue << " and gloablMaxValue is " << gloablMaxValue << " estimated err is " << err * (gloablMaxValue - gloablMinValue) << std::endl;
    vtkm::Float64 updatedError = err * (gloablMaxValue - gloablMinValue);
    // using pointNeighborhood worklet to process the data
    timer.Start();
    callCriticalPointWorklet(inData, updatedError, numSamples);
    timer.Stop();
    std::cout << "filter execution time: " << timer.GetElapsedTime() << std::endl;

    std::string outputFileName = "MinProb_MC_Weather" + std::to_string(dimx) + "_" + std::to_string(dimy) + "_" + std::to_string(numSamples) + ".vtk";
    vtkm::io::VTKDataSetWriter writeCross(outputFileName);
    writeCross.WriteDataSet(inData);

    return 0;
}