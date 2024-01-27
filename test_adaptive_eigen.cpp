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

#include "ucvworklet/ExtractingMeanStdev.hpp"
#include "ucvworklet/EntropyAdaptiveEigens.hpp"

#include <vtkm/cont/Timer.h>

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

void ComputeEntropyWithRuntimeVec(vtkm::cont::DataSet vtkmDataSet,
                                  double isovalue, int numSamples, std::string outputFileNameSuffix, bool use2d)
{
    // Processing current ensemble data sets based on uncertianty countour
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
    auto resolveType = [&](auto &concreteArray)
    {
        vtkm::cont::Invoker invoke;
        vtkm::Id numPoints = concreteArray.GetNumberOfValues();
        auto concreteArrayView = vtkm::cont::make_ArrayHandleView(concreteArray, 0, numPoints);

        invoke(ExtractingMean{}, concreteArrayView, meanArray);
        // printSummary_ArrayHandle(meanArray, std::cout);
        // printSummary_ArrayHandle(stdevArray, std::cout);

        vtkm::cont::ArrayHandle<vtkm::Float64> crossProbability;
        vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
        vtkm::cont::ArrayHandle<vtkm::Float64> entropy;

        // check 2d or 3d
        if (use2d)
        {
            invoke(EntropyAdaptiveEigens<4, 16>{isovalue, numSamples, 1.0, 0.0}, vtkmDataSet.GetCellSet(), concreteArrayView, meanArray, crossProbability, numNonZeroProb, entropy);
        }
        else
        {
            invoke(EntropyAdaptiveEigens<8, 256>{isovalue, numSamples, 1.0, 0.1}, vtkmDataSet.GetCellSet(), concreteArrayView, meanArray, crossProbability, numNonZeroProb, entropy);
        }

        auto outputDataSet = vtkmDataSet;

        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << isovalue;
        std::string isostr = stream.str();

        std::string outputFileName = outputFileNameSuffix + isostr + ".vtk";

        outputDataSet.AddCellField("cross_prob_" + isostr, crossProbability);
        outputDataSet.AddCellField("num_nonzero_prob" + isostr, numNonZeroProb);
        outputDataSet.AddCellField("entropy" + isostr, entropy);

        vtkm::io::VTKDataSetWriter writeCross(outputFileName);
        writeCross.WriteDataSet(outputDataSet);
    };

    vtkmDataSet.GetField("ensembles")
        .GetData()
        .CastAndCallWithExtractedArray(resolveType);
}

int main(int argc, char *argv[])
{
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    vtkm::cont::Timer timer{initResult.Device};

    if (argc != 10)
    {
        //./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
        std::cout << "<executable> <SyntheticDataSuffix> <FieldName> <Dimx> <Dimy> <Dimz> <iso> <num of sampls for mv> <num of ensembles> <outputFileSuffix>" << std::endl;
        exit(0);
    }

    std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

    std::string dataPathSuffix = std::string(argv[1]);
    std::string fieldName = std::string(argv[2]);

    int dimx = std::stoi(argv[3]);
    int dimy = std::stoi(argv[4]);
    int dimz = std::stoi(argv[5]);

    double isovalue = std::stod(argv[6]);
    int numSamples = std::stoi(argv[7]);
    int totalNumEnsemble = std::stoi(argv[8]);
    std::string outputSuffix = std::string(argv[9]);

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

    std::cout << "ok to load the data at the first step" << std::endl;

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
    bool use2d = false;
    if (dimz == 1)
    {
        use2d = true;
    }
    ComputeEntropyWithRuntimeVec(vtkmDataSet, isovalue, numSamples, outputSuffix, use2d);
    std::cout << "ok to get entropy for all ensembles" << std::endl;

    return 0;
}
