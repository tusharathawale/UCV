#include <float.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/DispatcherPointNeighborhood.h>
#include "ucvworklet/ExtractingByNeigoborhoodMinMax.hpp"
#include "ucvworklet/EntropyUniform.hpp"

#include "ucvworklet/ExtractingByNeigoborhoodMeanStdev.hpp"
#include "ucvworklet/EntropyIndependentGaussian.hpp"

#include "ucvworklet/ExtractingByNeigoborhoodMeanRaw.hpp"
#include "ucvworklet/MVGaussianWithEnsemble3DTryLialg.hpp"

#include <sstream>
#include <iomanip>

int oneDBlocks = 16;
int threadsPerBlock = 16;
#ifdef VTKM_CUDA
vtkm::cont::cuda::ScheduleParameters
mySchedParams(char const *name,
              int major,
              int minor,
              int multiProcessorCount,
              int maxThreadsPerMultiProcessor,
              int maxThreadsPerBlock)
{
    vtkm::cont::cuda::ScheduleParameters p;
    p.one_d_blocks = oneDBlocks;
    p.one_d_threads_per_block = threadsPerBlock;

    return p;
}
#endif

int main(int argc, char *argv[])
{

    // init the vtkm (set the backend and log level here)
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    vtkm::cont::Timer timer{initResult.Device};
    vtkm::cont::Timer timerDetails{initResult.Device};

    std::cout << "initResult.Device: " << initResult.Device.GetName() << " timer device: " << timer.GetDevice().GetName() << std::endl;

    if (argc != 7)
    {
        std::cout << "executable [VTK-m options] <filename> <fieldname> <distribution> <blocksize> <isovalue> <numSamples>" << std::endl;
        std::cout << "VTK-m options are:\n";
        std::cout << initResult.Usage << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    std::string distribution = argv[3];
    int blocksize = std::stoi(argv[4]);
    double isovalue = std::atof(argv[5]);
    int numSamples = std::atof(argv[6]);

    // this might cause some kokkos finalize error
    // since this global vtkm object is created before the init operation
    // of the kokkos
    using SupportedTypes = vtkm::List<vtkm::Float32,
                                      vtkm::Float64,
                                      vtkm::Int8,
                                      vtkm::UInt8,
                                      vtkm::Int16,
                                      vtkm::UInt16,
                                      vtkm::Int32,
                                      vtkm::UInt32,
                                      vtkm::Id>;

    /*

    #ifdef VTKM_CUDA

        if (backend == "cuda")
        {
            char const *nblock = getenv("UCV_GPU_NUMBLOCK");
            char const *nthread = getenv("UCV_GPU_BLOCKPERTHREAD");
            if (nblock != NULL && nthread != NULL)
            {
                oneDBlocks = std::stoi(std::string(nblock));
                threadsPerBlock = std::stoi(std::string(nthread));
                // the input value for the init scheduled parameter is a function
                vtkm::cont::cuda::InitScheduleParameters(mySchedParams);
                std::cout << "cuda parameters: " << oneDBlocks << " " << threadsPerBlock << std::endl;
            }
        }

    #endif
    */
    // load the dataset (beetles data set, structured one)
    // TODO, the data set can be distributed between different ranks

    // create the vtkm data set from the loaded data
    std::cout << "fileName: " << fileName << std::endl;
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // check the property of the data
    inData.PrintSummary(std::cout);

    // TODO timer start to extract key
    // auto timer1 = std::chrono::steady_clock::now();

    timer.Synchronize();
    timer.Start();
    timerDetails.Start();

    auto field = inData.GetField(fieldName);

    auto cellSet = inData.GetCellSet();

    // Assuming the imput data is the structured data

    bool isStructured = cellSet.IsType<vtkm::cont::CellSetStructured<3>>();
    if (!isStructured)
    {
        std::cout << "the extraction only works for CellSetStructured<3>" << std::endl;
        exit(0);
    }

    vtkm::cont::CellSetStructured<3> structCellSet =
        cellSet.AsCellSet<vtkm::cont::CellSetStructured<3>>();

    vtkm::Id3 pointDims = structCellSet.GetPointDimensions();

    std::cout << "------" << std::endl;
    std::cout << "point dim: " << pointDims[0] << " " << pointDims[1] << " " << pointDims[2] << std::endl;

    // go through all points and set the specific key
    vtkm::Id xdim = pointDims[0];
    vtkm::Id ydim = pointDims[1];
    vtkm::Id zdim = pointDims[2];

    if (xdim % blocksize != 0 || ydim % blocksize != 0 || zdim % blocksize != 0)
    {
        throw std::runtime_error("dim is supposed to be dividied by blocksize");
    }

    vtkm::Id numberBlockx = xdim % blocksize == 0 ? xdim / blocksize : xdim / blocksize + 1;
    vtkm::Id numberBlocky = ydim % blocksize == 0 ? ydim / blocksize : ydim / blocksize + 1;
    vtkm::Id numberBlockz = zdim % blocksize == 0 ? zdim / blocksize : zdim / blocksize + 1;

    const vtkm::Id3 reducedDims(numberBlockx, numberBlocky, numberBlockz);

    std::cout << "reducedDims " << reducedDims << std::endl;
    timerDetails.Stop();
    std::cout << "timerDetails 1 " << timerDetails.GetElapsedTime() * 1000 << std::endl;
    timerDetails.Start();

    //auto coords = inData.GetCoordinateSystem();
    auto bounds = inData.GetCoordinateSystem().GetBounds();
    timerDetails.Stop();
    std::cout << "timerDetails 1.1 " << timerDetails.GetElapsedTime() * 1000 << std::endl;
    timerDetails.Start();
    
    auto reducedOrigin = bounds.MinCorner();

    vtkm::FloatDefault spacex = (bounds.X.Max - bounds.X.Min) / (numberBlockx - 1);
    vtkm::FloatDefault spacey = (bounds.Y.Max - bounds.Y.Min) / (numberBlocky - 1);
    vtkm::FloatDefault spacez = (bounds.Z.Max - bounds.Z.Min) / (numberBlockz - 1);

    vtkm::Vec3f_64 reducedSpaceing(spacex, spacey, spacez);

    timerDetails.Stop();
    std::cout << "timerDetails 1.2 " << timerDetails.GetElapsedTime() * 1000 << std::endl;
    timerDetails.Start();


    timerDetails.Stop();
    std::cout << "timerDetails 1.3 " << timerDetails.GetElapsedTime() * 1000 << std::endl;
    timerDetails.Start();

    vtkm::cont::DataSetBuilderUniform dataSetBuilder;
    // origin is {0,0,0} spacing is {blocksize,blocksize,blocksize} make sure the reduced data
    // are in same shape with original data
    vtkm::cont::DataSet reducedDataSet = dataSetBuilder.Create(reducedDims, reducedOrigin, reducedSpaceing);

    timerDetails.Stop();
    std::cout << "timerDetails 2 " << timerDetails.GetElapsedTime() * 1000 << std::endl;
    timerDetails.Start();

    // declare results array
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> entropyResult;

    if (distribution == "uni")
    {
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> minArray;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxArray;

        // using the cellset of reducedDataSet
        // using the field data from original data
        // only one pass
        using WorkletTypeNMinMax = ExtractingByNeigoborhoodMinMax;
        using DispatcherType = vtkm::worklet::DispatcherPointNeighborhood<WorkletTypeNMinMax>;

        auto resolveType = [&](const auto &concrete)
        {
            DispatcherType dispatcher(WorkletTypeNMinMax{isovalue, numSamples, blocksize, xdim, ydim, zdim});
            dispatcher.Invoke(reducedDataSet.GetCellSet(), concrete, minArray, maxArray);
        };

        field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
            resolveType);
        timer.Synchronize();
        timer.Stop();

        timerDetails.Stop();
        std::cout << "timerDetails 3 " << timerDetails.GetElapsedTime() * 1000 << std::endl;

        std::cout << "sampling min and max time " << timer.GetElapsedTime() * 1000 << std::endl;

        timer.Synchronize();
        timer.Start();

        using WorkletType = EntropyUniform;
        using DispatcherEntropyUniform = vtkm::worklet::DispatcherMapTopology<WorkletType>;

        DispatcherEntropyUniform dispatcherEntropyUniform(EntropyUniform{isovalue});
        dispatcherEntropyUniform.Invoke(reducedDataSet.GetCellSet(), minArray, maxArray, crossProb, numNonZeroProb, entropyResult);
        timer.Synchronize();
        timer.Stop();

        std::cout << "uncertainty uni time " << timer.GetElapsedTime() * 1000 << std::endl;
    }
    else if (distribution == "ig")
    {
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdevArray;

        using WorkletTypeNMeanStdev = ExtractingByNeigoborhoodMeanStdev;
        using DispatcherType = vtkm::worklet::DispatcherPointNeighborhood<WorkletTypeNMeanStdev>;

        auto resolveType = [&](const auto &concrete)
        {
            DispatcherType dispatcher(WorkletTypeNMeanStdev{isovalue, numSamples, blocksize, xdim, ydim, zdim});
            dispatcher.Invoke(reducedDataSet.GetCellSet(), concrete, meanArray, stdevArray);
        };

        field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
            resolveType);
        timer.Synchronize();
        timer.Stop();

        std::cout << "sampling mean and stdev time " << timer.GetElapsedTime() * 1000 << std::endl;
        timer.Synchronize();
        timer.Start();
        using WorkletType = EntropyIndependentGaussian;
        using DispatcherEntropyIG = vtkm::worklet::DispatcherMapTopology<WorkletType>;

        DispatcherEntropyIG dispatcherEntropyIG(EntropyIndependentGaussian{isovalue});
        dispatcherEntropyIG.Invoke(reducedDataSet.GetCellSet(), meanArray, stdevArray, crossProb, numNonZeroProb, entropyResult);
        timer.Synchronize();
        timer.Stop();
        std::cout << "EIGaussianTime time: " << timer.GetElapsedTime() * 1000 << std::endl;
    }
    else if (distribution == "mg")
    {
        if (xdim % 4 != 0 || ydim % 4 != 0 || zdim % 4 != 0)
        {
            // if the data size is not divided by blocksize
            // we can reample or padding the data set before hand
            // it will be convenient to compute cov matrix by this way
            throw std::runtime_error("only support blocksize = 4 and the case where xyz dim is diveide dy blocksize for current mg");
        }

        using VecType = vtkm::Vec<vtkm::FloatDefault, 4 * 4 * 4>;
        vtkm::cont::ArrayHandle<VecType> SOARawArray;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;

        using WorkletTypeMeanRaw = ExtractingByNeigoborhoodMeanRaw;
        using DispatcherType = vtkm::worklet::DispatcherPointNeighborhood<WorkletTypeMeanRaw>;

        auto resolveType = [&](const auto &concrete)
        {
            DispatcherType dispatcher(WorkletTypeMeanRaw{isovalue, numSamples, blocksize, xdim, ydim, zdim});
            dispatcher.Invoke(reducedDataSet.GetCellSet(), concrete, meanArray, SOARawArray);
        };

        field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
            resolveType);
        timer.Synchronize();
        timer.Stop();

        std::cout << "sampling mean and raw time " << timer.GetElapsedTime() * 1000 << std::endl;
        timer.Synchronize();
        timer.Start();

        using WorkletTypeMVG = MVGaussianWithEnsemble3DTryLialg;
        using DispatcherTypeMVG = vtkm::worklet::DispatcherMapTopology<WorkletTypeMVG>;

        DispatcherTypeMVG dispatcherMVG(WorkletTypeMVG{isovalue, numSamples});
        dispatcherMVG.Invoke(reducedDataSet.GetCellSet(), SOARawArray, meanArray, crossProb, numNonZeroProb, entropyResult);
        timer.Synchronize();
        timer.Stop();
        std::cout << "MVGTime time: " << timer.GetElapsedTime() * 1000 << std::endl;
    }
    else
    {
        throw std::runtime_error("unsupported distribution: " + distribution);
    }

    /* write results*/
    reducedDataSet.AddCellField("entropy", entropyResult);
    reducedDataSet.AddCellField("num_nonzero_prob", numNonZeroProb);
    reducedDataSet.AddCellField("cross_prob", crossProb);

    // using the same type as the assumption for the output type
    // std::cout << "===data summary for reduced data with uncertainty:" << std::endl;

    // reducedDataSet.PrintSummary(std::cout);
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << isovalue;
    std::string isostr = stream.str();

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + "_iso" + isostr + "_" + distribution + "_block" + std::to_string(blocksize) + std::string("_NeighborProb.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(reducedDataSet);

    return 0;
}
