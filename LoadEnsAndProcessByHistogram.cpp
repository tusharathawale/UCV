#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>

#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandleRuntimeVec.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <CriticalPointsHistogram.h>
#include <EnsembleFieldToHistogramField.h>

int main(int argc, char *argv[])
{
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    vtkm::cont::Timer timer{initResult.Device};

    if (argc != 8)
    {
        //./test_syntheticdata_el_sequence /Users/zw1/Documents/cworkspace/src/UCV/exp_scripts/create_dataset/RawdataPointScalar TestField 300 0.8 1000
        std::cout << "<executable> <SyntheticDataSuffix> <FieldName> <Dimx> <Dimy> <Dimz> <num of ensembles> <num of bins>" << std::endl;
        exit(0);
    }

    std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;

    std::string dataPathSuffix = std::string(argv[1]);
    std::string fieldName = std::string(argv[2]);

    int dimx = std::stoi(argv[3]);
    int dimy = std::stoi(argv[4]);
    int dimz = std::stoi(argv[5]);

    int numEnsembles = std::stoi(argv[6]);
    int numHistBin = std::stoi(argv[7]);

    const vtkm::Id3 dims(dimx, dimy, dimz);
    vtkm::cont::DataSetBuilderUniform dataSetBuilder;
    vtkm::cont::DataSet ensembleData = dataSetBuilder.Create(dims);

    // load all ens data and put them into one entry
    // for the large ens array
    // the length of runtime vec is dimx * dimy * dimz
    // each element is an vector with numEnsemble values
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> allEnsemblesArray(numEnsembles);
    allEnsemblesArray.Allocate(dimx * dimy * dimz);

    std::vector<vtkm::cont::ArrayHandle<vtkm::FloatDefault>> dataArray;
    // redsea data start from 1
    for (int ensId = 0; ensId < numEnsembles; ensId++)
    {
        std::string fileName = dataPathSuffix + "_" + std::to_string(ensId) + ".vtk";
        vtkm::io::VTKDataSetReader reader(fileName);
        vtkm::cont::DataSet inData = reader.ReadDataSet();

        vtkm::cont::ArrayHandle<vtkm::FloatDefault> fieldDataArray;
        vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField(fieldName).GetData(), fieldDataArray);
        dataArray.push_back(fieldDataArray);
    }

    //std::cout << "ok to load the data at the first step" << std::endl;

    // using all ensembles
    vtkm::cont::ArrayHandleRuntimeVec<vtkm::FloatDefault> runtimeVecArray(numEnsembles);
    runtimeVecArray.Allocate(dimx * dimy);
    auto writePortal = runtimeVecArray.WritePortal();
    for (int j = 0; j < dimy; j++)
    {
        for (int i = 0; i < dimx; i++)
        {
            int pointIndex = j * dimx + i;
            auto vecValue = writePortal.Get(pointIndex);
            // load ensemble data from ens 0 to ens with id usedEnsembles-1
            for (int currEndId = 0; currEndId < numEnsembles; currEndId++)
            {
                // set ensemble value
                vecValue[currEndId] = dataArray[currEndId].ReadPortal().Get(pointIndex);
            }
        }
    }

    ensembleData.AddPointField("ensembles", runtimeVecArray);
    // printSummary_ArrayHandle(runtimeVecArray, std::cout);

    timer.Start();
    // Should we be timing the conversion of ensemble data to uncertain data?
    vtkm::filter::uncertainty::EnsembleFieldToHistogramField ensemble2histogram;
    ensemble2histogram.SetNumberOfBins(numHistBin); // Right amount?
    ensemble2histogram.SetActiveField("ensembles", vtkm::cont::Field::Association::Points);
    ensemble2histogram.SetHistogramDensityName("HistogramDensity");
    ensemble2histogram.SetHistogramEdgesName("HistogramEdges");
    vtkm::cont::DataSet uncertainData = ensemble2histogram.Execute(ensembleData);

    vtkm::filter::uncertainty::CriticalPointsHistogram criticalPoints;
    criticalPoints.SetNumberOfBins(numHistBin);
    criticalPoints.SetHistogramDensityName("HistogramDensity");
    criticalPoints.SetHistogramEdgesName("HistogramEdges");
    criticalPoints.SetMinimumProbabilityName("MinimumProbability");
    vtkm::cont::DataSet criticalPointsData = criticalPoints.Execute(uncertainData);
    timer.Stop();
    std::cout << "execution time of critical points filter is " << timer.GetElapsedTime() << std::endl;

    std::string outputFileName = "MinProb_hist_" + std::to_string(dimx) + "_" + std::to_string(dimy) + ".vtk";
    vtkm::io::VTKDataSetWriter writeCross(outputFileName);
    writeCross.WriteDataSet(criticalPointsData);

    return 0;
}
