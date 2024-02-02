#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandleRuntimeVec.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/Algorithm.h>

#include "ucvworklet/ExtractingMeanStdev.hpp"
#include "ucvworklet/EntropyAdaptiveEigensPearson.hpp"
#include "ucvworklet/EntropyAdaptiveEigens.hpp"

#include "ucvworklet/MVGaussianWithEnsemble3DTryEL.hpp"
#include "ucvworklet/ComputeDiffSum.hpp"

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

// compute the RMSE

vtkm::cont::ArrayHandle<vtkm::Float64> ComputeEntropyWithRuntimeVec(vtkm::cont::DataSet vtkmDataSet,
                                                                    double isovalue, int numSamples, std::string outputFileNameSuffix, bool use2d, double eigenThreshold, vtkm::cont::Timer &timer, std::string writeFile)
{
    timer.Start();
    // Processing current ensemble data sets based on uncertianty countour
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdArray;

    vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::Float64> entropy;
    auto resolveType = [&](auto &concreteArray)
    {
        vtkm::cont::Invoker invoke;
        vtkm::Id numPoints = concreteArray.GetNumberOfValues();
        auto concreteArrayView = vtkm::cont::make_ArrayHandleView(concreteArray, 0, numPoints);

        invoke(ExtractingMeanStdevEnsembles{}, concreteArrayView, meanArray, stdArray);
        // printSummary_ArrayHandle(meanArray, std::cout);
        // printSummary_ArrayHandle(stdevArray, std::cout);

        // check 2d or 3d
        if (use2d)
        {
            invoke(EntropyAdaptiveEigens<4, 16>{isovalue, numSamples, 1.0, eigenThreshold}, vtkmDataSet.GetCellSet(), concreteArrayView, meanArray, crossProbability, numNonZeroProb, entropy);
        }
        else
        {
            invoke(EntropyAdaptiveEigens<8, 256>{isovalue, numSamples, 1.0, eigenThreshold}, vtkmDataSet.GetCellSet(), concreteArrayView, meanArray, crossProbability, numNonZeroProb, entropy);
        }

        auto outputDataSet = vtkmDataSet;

        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << isovalue;
        std::string isostr = stream.str();

        std::string outputFileName = outputFileNameSuffix + isostr + ".vtk";

        outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);
        outputDataSet.AddCellField("num_nonzero_prob" + isostr, numNonZeroProb);
        outputDataSet.AddCellField("entropy" + isostr, entropy);
        if (writeFile == "true")
        {
            vtkm::io::VTKDataSetWriter writeCross(outputFileName);
            writeCross.WriteDataSet(outputDataSet);
        }
    };

    vtkmDataSet.GetField("ensembles")
        .GetData()
        .CastAndCallWithExtractedArray(resolveType);
    timer.Stop();
    std::cout << "worklet time is:" << timer.GetElapsedTime() << std::endl;

    return crossProbability;
}

vtkm::cont::ArrayHandle<vtkm::Float64> ComputeEntropyBasedOnPearson(vtkm::cont::DataSet vtkmDataSet, vtkm::cont::DataSet pearsonData,
                                                                    double isovalue, int numSamples, std::string outputFileNameSuffix, double eigenThreshold, double pearsonThreshold,vtkm::cont::Timer &timer)
{
    timer.Start();
    // Processing current ensemble data sets based on uncertianty countour
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
    // add std
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdArray;

    vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::Float64> entropy;
    vtkm::cont::ArrayHandle<vtkm::Float64> pearsonCorrelation;
    //pearsonData.PrintSummary(std::cout);
    vtkm::cont::ArrayCopyShallowIfPossible(pearsonData.GetField("maxPearsonCorrelation").GetData(), pearsonCorrelation);
    auto resolveType = [&](auto &concreteArray)
    {
        vtkm::cont::Invoker invoke;
        vtkm::Id numPoints = concreteArray.GetNumberOfValues();
        auto concreteArrayView = vtkm::cont::make_ArrayHandleView(concreteArray, 0, numPoints);
        invoke(ExtractingMeanStdevEnsembles{}, concreteArrayView, meanArray, stdArray);

        //std::cout << "pearson correlation values: " << pearsonCorrelation.GetNumberOfValues() << std::endl;
        //std::cout << "ensmeble data set: " << vtkmDataSet.GetNumberOfCells() << std::endl;
        if(pearsonCorrelation.GetNumberOfValues()!=vtkmDataSet.GetNumberOfCells()){
            throw std::runtime_error("pearsonCorrelation values is not equal to vtkmDataSet cells");
        }

        //visit cell based on pearson correlations
        invoke(EntropyAdaptiveEigensPearson<8, 256>{isovalue, numSamples, 1.0, eigenThreshold, pearsonThreshold}, vtkmDataSet.GetCellSet(), concreteArray, meanArray, stdArray, pearsonCorrelation, crossProbability, numNonZeroProb, entropy);


    };

    vtkmDataSet.GetField("ensembles")
        .GetData()
        .CastAndCallWithExtractedArray(resolveType);
    timer.Stop();
    std::cout << "worklet time is:" << timer.GetElapsedTime() << std::endl;

    return crossProbability;
}

vtkm::cont::ArrayHandle<vtkm::Float64> ComputeEntropyWithOrigianlMVG(vtkm::cont::DataSet vtkmDataSet,
                                                                     double isovalue, int numSamples, std::string outputFileNameSuffix, bool use2d, double eigenThreshold, vtkm::cont::Timer &timer)
{
    timer.Start();
    // Processing current ensemble data sets based on uncertianty countour
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
    vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::Float64> entropy;
    auto resolveType = [&](auto &concreteArray)
    {
        vtkm::cont::Invoker invoke;
        vtkm::Id numPoints = concreteArray.GetNumberOfValues();
        auto concreteArrayView = vtkm::cont::make_ArrayHandleView(concreteArray, 0, numPoints);

        invoke(ExtractingMean{}, concreteArrayView, meanArray);
        // printSummary_ArrayHandle(meanArray, std::cout);
        // printSummary_ArrayHandle(stdevArray, std::cout);

        // check 2d or 3d
        if (use2d)
        {
            std::runtime_error("only support 3d case for this testing");
        }
        else
        {
            invoke(MVGaussianWithEnsemble3DTryEL{isovalue, numSamples}, vtkmDataSet.GetCellSet(), concreteArrayView, meanArray, crossProbability, numNonZeroProb, entropy);
        }

        auto outputDataSet = vtkmDataSet;

        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << isovalue;
        std::string isostr = stream.str();

        std::string outputFileName = outputFileNameSuffix + isostr + ".vtk";

        outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);
        outputDataSet.AddCellField("num_nonzero_prob" + isostr, numNonZeroProb);
        outputDataSet.AddCellField("entropy" + isostr, entropy);
    };

    vtkmDataSet.GetField("ensembles")
        .GetData()
        .CastAndCallWithExtractedArray(resolveType);
    timer.Stop();
    std::cout << "worklet time is:" << timer.GetElapsedTime() << std::endl;

    return crossProbability;
}

int main(int argc, char *argv[])
{
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    vtkm::cont::Timer timer{initResult.Device};

    if (argc != 13)
    {
        // ./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
        std::cout << "<executable> <SyntheticDataSuffix> <PearsonFile> <FieldName> <Dimx> <Dimy> <Dimz> <iso> <num of sampls for mv> <num of ensembles> <outputFileSuffix> <eigenThreshold> <pearsonThreshold>" << std::endl;
        exit(0);
    }

    std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

    std::string dataPathSuffix = std::string(argv[1]);
    std::string PearsonFile = std::string(argv[2]);
    std::string fieldName = std::string(argv[3]);

    int dimx = std::stoi(argv[4]);
    int dimy = std::stoi(argv[5]);
    int dimz = std::stoi(argv[6]);

    double isovalue = std::stod(argv[7]);
    int numSamples = std::stoi(argv[8]);
    int totalNumEnsemble = std::stoi(argv[9]);
    std::string outputSuffix = std::string(argv[10]);
    double eigenThreshold = std::stod(argv[11]);
    double pearsonThreshold = std::stod(argv[12]);

    if (eigenThreshold < 0)
    {
        throw std::runtime_error("eigenThreshold is supposed to be larger than 0");
    }

    std::cout << "iso value is: " << isovalue << " numSamples is: " << numSamples << std::endl;

    vtkm::Id xdim = dimx;
    vtkm::Id ydim = dimy;
    vtkm::Id zdim = dimz;

    const vtkm::Id3 dims(xdim, ydim, zdim);
    vtkm::cont::DataSetBuilderUniform dataSetBuilder;
    vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);

    // load all data values
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> allEnsemblesArray(totalNumEnsemble);
    allEnsemblesArray.Allocate(dimx * dimy * dimz);

    // for each ensemble(version) of the data
    // store results into the allEnsemblesArray
    std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> dataArray;
    for (int ensId = 0; ensId < totalNumEnsemble; ensId++)
    {
        std::string fileName = dataPathSuffix + "_" + std::to_string(ensId) + ".vtk";
        vtkm::io::VTKDataSetReader reader(fileName);
        vtkm::cont::DataSet inData = reader.ReadDataSet();

        vtkm::cont::ArrayHandle<vtkm::Float64> fieldDataArray;
        vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(fieldName).GetData(), fieldDataArray);
        dataArray.push_back(fieldDataArray);
    }

    std::cout << "ok to load the data" << std::endl;

    // using all ensembles
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> runtimeVecArray(totalNumEnsemble);
    runtimeVecArray.Allocate(dimx * dimy * dimz);
    auto writePortal = runtimeVecArray.WritePortal();
    for (int k = 0; k < dimz; k++)
    {
        for (int j = 0; j < dimy; j++)
        {
            for (int i = 0; i < dimx; i++)
            {
                int pointIndex = k * dimx * dimy + j * dimx + i;
                auto vecValue = writePortal.Get(pointIndex);
                // load ensemble data from ens 0 to ens with id usedEnsembles-1
                for (int currEndId = 0; currEndId < totalNumEnsemble; currEndId++)
                {
                    // set ensemble value
                    vecValue[currEndId] = dataArray[currEndId].ReadPortal().Get(pointIndex);
                }
            }
        }
    }

    vtkmDataSet.AddPointField("ensembles", runtimeVecArray);

    // load the pearson dataset
    vtkm::io::VTKDataSetReader pearsonReader(PearsonFile);
    vtkm::cont::DataSet pearsonData = pearsonReader.ReadDataSet();

    std::cout << "start to call the data by pearson correlation" << std::endl;
    auto crossProb1 = ComputeEntropyBasedOnPearson(vtkmDataSet, pearsonData, isovalue, numSamples, outputSuffix + "_pearson", eigenThreshold, pearsonThreshold, timer);

    std::cout << "start to call the original mvg worklet" << std::endl;
    auto crossProb2 = ComputeEntropyWithOrigianlMVG(vtkmDataSet, isovalue, numSamples, outputSuffix + "_originalMVG", false, eigenThreshold, timer);

    // compute the RMSE for two cases
    vtkm::cont::Invoker invoke;
    vtkm::cont::ArrayHandle<vtkm::Float64> diff;

    invoke(ComputeDiffSquare{}, crossProb1, crossProb2, diff);
    vtkm::Float64 diffSquarSum = vtkm::cont::Algorithm::Reduce(diff, 0.0, vtkm::Sum());
    
    std::cout << "DiffSquarSum is: " << diffSquarSum << " RMSE is: " << vtkm::Sqrt(diffSquarSum/diff.GetNumberOfValues()) << std::endl;


    return 0;
}
