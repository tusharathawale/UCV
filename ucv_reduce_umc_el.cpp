#include <float.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/DispatcherReduceByKey.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/DispatcherMapTopology.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/ExtractingMinMax.hpp"
#include "ucvworklet/ExtractingMeanStdev.hpp"

#include "ucvworklet/EntropyUniform.hpp"
#include "ucvworklet/EntropyIndependentGaussian.hpp"

#include <sstream>
#include <iomanip>

// #include <chrono>

// #ifdef VTKM_CUDA
// #else
//  the case three does not works well for cuda at this step
#include "ucvworklet/ExtractingMeanRaw.hpp"

#include "ucvworklet/MVGaussianWithEnsemble3DTryEL.hpp"

// #endif // VTKM_CUDA

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

std::string backend = "openmp";

void initBackend(vtkm::cont::Timer &timer)
{
    // init the vtkh device
    char const *tmp = getenv("UCV_VTKM_BACKEND");

    if (tmp == nullptr)
    {
        return;
    }
    else
    {
        backend = std::string(tmp);
        std::cout << "Setting the device with UCV_VTKM_BACKEND=" << backend << "\n";
        std::cout << "This method is antiquated. Consider using the --vtkm-device command line argument." << std::endl;
    }

    // if (rank == 0)
    //{
    std::cout << "vtkm backend is:" << backend << std::endl;
    //}

    if (backend == "serial")
    {
        vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagSerial());
        timer.Reset(vtkm::cont::DeviceAdapterTagSerial());
    }
    else if (backend == "openmp")
    {
        vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP());
        timer.Reset(vtkm::cont::DeviceAdapterTagOpenMP());
    }
    else if (backend == "cuda")
    {
        vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagCuda());
        timer.Reset(vtkm::cont::DeviceAdapterTagCuda());
    }
    else
    {
        std::cerr << " unrecognized backend " << backend << std::endl;
    }
    return;
}

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

int main(int argc, char *argv[])
{

    // init the vtkm (set the backend and log level here)
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    vtkm::cont::Timer timer{initResult.Device};
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

    auto keyArray =
        vtkm::cont::ArrayHandleCounting<vtkm::Id>(0, 1, static_cast<vtkm::Id>(xdim * ydim * zdim));

    vtkm::Id numberBlockx = xdim % blocksize == 0 ? xdim / blocksize : xdim / blocksize + 1;
    vtkm::Id numberBlocky = ydim % blocksize == 0 ? ydim / blocksize : ydim / blocksize + 1;
    vtkm::Id numberBlockz = zdim % blocksize == 0 ? zdim / blocksize : zdim / blocksize + 1;

    // add key array into the dataset, and check the output
    // inData.AddPointField("keyArray", keyArrayNew);
    // std::cout << "------" << std::endl;
    // inData.PrintSummary(std::cout);

    // std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    // std::string outputFileName = fileSuffix + std::string("_Key.vtk");
    // vtkm::io::VTKDataSetWriter write(outputFileName);
    // write.WriteDataSet(inData);

    // TODO, the decision of the distribution should start from this position
    // for uniform case, we extract min and max
    // for gaussian case, we extract other values

    // create the new data sets for the reduced data
    // the dims for new data sets are numberBlockx*numberBlocky*numberBlockz
    const vtkm::Id3 reducedDims(numberBlockx, numberBlocky, numberBlockz);

    // vtkm introduces some redoundant operations here
    // https://gitlab.kitware.com/vtk/vtk-m/-/merge_requests/3008/pipelines
    // we just comment out this part temporarily and hard code the value
    //  to make sure the speed up number makes sense
    //  auto bounds = inData.GetCoordinateSystem().GetBounds();

    // this only works for the beetle data, for the temporary evaluation
    // we need to add things back for other data sets
    // { X:[0..495], Y:[0..831], Z:[0..831] }
    // another way is to
    // get the coordinates as an `ArrayHandleUniformPointCoordinates`
    // and get the bounds from that
    // this can also avoid the bounds issue
    // beetle
    vtkm::Bounds bounds(0,495,0,831,0,831);
    // supernova
    // vtkm::Bounds bounds(-1.02655e+09, 9.42128e+08, -1.02353e+09, 1.00073e+09, -9.16584e+08, 9.32988e+08);

    auto reducedOrigin = bounds.MinCorner();

    vtkm::FloatDefault spacex = (bounds.X.Max - bounds.X.Min) / (numberBlockx - 1);
    vtkm::FloatDefault spacey = (bounds.Y.Max - bounds.Y.Min) / (numberBlocky - 1);
    vtkm::FloatDefault spacez = (bounds.Z.Max - bounds.Z.Min) / (numberBlockz - 1);

    vtkm::Vec3f_64 reducedSpaceing(spacex, spacey, spacez);

    vtkm::cont::DataSetBuilderUniform dataSetBuilder;
    // origin is {0,0,0} spacing is {blocksize,blocksize,blocksize} make sure the reduced data
    // are in same shape with original data
    vtkm::cont::DataSet reducedDataSet = dataSetBuilder.Create(reducedDims, reducedOrigin, reducedSpaceing);

    // declare results array
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProb;
    vtkm::cont::ArrayHandle<vtkm::Id> numNonZeroProb;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> entropyResult;

    // Step1 creating new key
    vtkm::cont::ArrayHandle<vtkm::Id> keyArrayNew;

    using DispatcherCreateKey = vtkm::worklet::DispatcherMapField<CreateNewKeyWorklet>;
    DispatcherCreateKey dispatcher(CreateNewKeyWorklet{xdim, ydim, zdim,
                                                       numberBlockx, numberBlocky, numberBlockz,
                                                       blocksize});

    dispatcher.Invoke(keyArray, keyArrayNew);

    // auto timer2 = std::chrono::steady_clock::now();
    // float extractKey =
    //     std::chrono::duration<float, std::milli>(timer2 - timer1).count();
    timer.Synchronize();
    timer.Stop();
    std::cout << "extractKey time: " << timer.GetElapsedTime() * 1000 << std::endl;
    timer.Synchronize();
    timer.Start();

    if (distribution == "uni")
    {

        // Step2 extracting ensemble data based on new key
        using DispatcherType = vtkm::worklet::DispatcherReduceByKey<ExtractingMinMax>;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> minArray;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxArray;
        vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

        auto resolveType = [&](const auto &concrete)
        {
            DispatcherType dispatcher;
            dispatcher.Invoke(keys, concrete, minArray, maxArray);
        };

        field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
            resolveType);

        // TODO timer ok to extract property
        // auto timer3 = std::chrono::steady_clock::now();
        // float extractMinMax =
        //    std::chrono::duration<float, std::milli>(timer3 - timer2).count();
        timer.Synchronize();
        timer.Stop();
        std::cout << "extractMinMax time: " << timer.GetElapsedTime() * 1000 << std::endl;
        timer.Synchronize();
        timer.Start();
        // generate the new data sets with min and max
        reducedDataSet.AddPointField("ensemble_min", minArray);
        reducedDataSet.AddPointField("ensemble_max", maxArray);
        reducedDataSet.PrintSummary(std::cout);

        // output the dataset into the vtk file for results checking
        // std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
        // std::string outputFileName = fileSuffix + "_" + fieldName + std::string("_ReduceMinMax.vtk");
        // vtkm::io::VTKDataSetWriter write(outputFileName);
        // write.WriteDataSet(reducedDataSet);

        // uniform distribution
        // using WorkletType = EntropyUniform;
        // using DispatcherEntropyUniform = vtkm::worklet::DispatcherMapTopology<WorkletType>;
        // DispatcherEntropyUniform dispatcherEntropyUniform(EntropyUniform{isovalue});
        // dispatcherEntropyUniform.Invoke(reducedDataSet.GetCellSet(), minArray, maxArray, crossProb, numNonZeroProb, entropyResult);

        // TODO make sure the worklet finish
        // TODO timer ok to compute uncertainty
        // auto timer4 = std::chrono::steady_clock::now();
        // float EntropyUniformTime =
        //    std::chrono::duration<float, std::milli>(timer4 - timer3).count();
        timer.Synchronize();
        timer.Stop();
        std::cout << "EntropyUniformTime time: " << timer.GetElapsedTime() * 1000 << std::endl;
    }
    else if (distribution == "ig")
    {
        // indepedent gaussian

        // extracting mean and stdev

        using DispatcherType = vtkm::worklet::DispatcherReduceByKey<ExtractingMeanStdev>;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> stdevArray;
        vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

        auto resolveType = [&](const auto &concrete)
        {
            DispatcherType dispatcher;
            dispatcher.Invoke(keys, concrete, meanArray, stdevArray);
        };

        field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
            resolveType);

        // TODO timer ok to extract properties
        // auto timer3 = std::chrono::steady_clock::now();
        // float ExtractingMeanStdevTime =
        //    std::chrono::duration<float, std::milli>(timer3 - timer2).count();
        timer.Synchronize();
        timer.Stop();
        std::cout << "ExtractingMeanStdev time: " << timer.GetElapsedTime() * 1000 << std::endl;
        timer.Synchronize();
        timer.Start();
        using WorkletType = EntropyIndependentGaussian<8,256>;
        using DispatcherEntropyIG = vtkm::worklet::DispatcherMapTopology<WorkletType>;

        DispatcherEntropyIG dispatcherEntropyIG(EntropyIndependentGaussian<8,256>{isovalue});
        dispatcherEntropyIG.Invoke(reducedDataSet.GetCellSet(), meanArray, stdevArray, crossProb, numNonZeroProb, entropyResult);

        // TODO timer ok to compute uncertianty
        // auto timer4 = std::chrono::steady_clock::now();
        // float EIGaussianTime =
        //    std::chrono::duration<float, std::milli>(timer4 - timer3).count();
        timer.Synchronize();
        timer.Stop();
        std::cout << "EIGaussianTime time: " << timer.GetElapsedTime() * 1000 << std::endl;
    }
    else if (distribution == "mg")
    {
        // #ifdef VTKM_CUDA
        //         std::cout << "multivariant gaussian does not work well for cuda now" << std::endl;
        //         exit(0);
        // #else
        //  multivariant gaussian
        //  extracting the mean and rawdata for each hixel block
        //  the raw data is used to compute the covariance matrix
        if (xdim % 4 != 0 || ydim % 4 != 0 || zdim % 4 != 0)
        {
            // if the data size is not divided by blocksize
            // we can reample or padding the data set before hand
            // it will be convenient to compute cov matrix by this way
            throw std::runtime_error("only support blocksize = 4 and the case where xyz dim is diveide dy blocksize for current mg");
        }
        // Step2 extracting the soa raw array
        // the value here should be same with the elements in each hixel

        using WorkletType = ExtractingMeanRaw;
        using DispatcherType = vtkm::worklet::DispatcherReduceByKey<WorkletType>;

        using VecType = vtkm::Vec<vtkm::FloatDefault, 4 * 4 * 4>;
        vtkm::cont::ArrayHandle<VecType> SOARawArray;
        vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
        // Pay attention to transfer the arrayHandle into the Keys type
        vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

        auto resolveType = [&](const auto &concrete)
        {
            DispatcherType dispatcher;
            dispatcher.Invoke(keys, concrete, meanArray, SOARawArray);
        };

        field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
            resolveType);

        // auto timer3 = std::chrono::steady_clock::now();
        // float ExtractingMeanRawTime =
        //     std::chrono::duration<float, std::milli>(timer3 - timer2).count();
        timer.Synchronize();
        timer.Stop();
        std::cout << "ExtractingMeanRawTime time: " << timer.GetElapsedTime() * 1000 << std::endl;
        timer.Synchronize();
        timer.Start();

        // step3 computing the cross probability
        // using WorkletTypeMVG = MVGaussianWithEnsemble3D;
        using WorkletTypeMVG = MVGaussianWithEnsemble3DTryEL;
        using DispatcherTypeMVG = vtkm::worklet::DispatcherMapTopology<WorkletTypeMVG>;

        DispatcherTypeMVG dispatcherMVG(WorkletTypeMVG{isovalue, numSamples});
        dispatcherMVG.Invoke(reducedDataSet.GetCellSet(), SOARawArray, meanArray, crossProb, numNonZeroProb, entropyResult);

        // TODO make sure it actually finish
        // auto timer4 = std::chrono::steady_clock::now();
        // float MVGTime =
        //    std::chrono::duration<float, std::milli>(timer4 - timer3).count();
        timer.Synchronize();
        timer.Stop();
        std::cout << "MVGTime time: " << timer.GetElapsedTime() * 1000 << std::endl;

        // #endif // VTKM_CUDA
    }
    else
    {
        throw std::runtime_error("unsupported distribution: " + distribution);
    }

    // using the same type as the assumption for the output type
    // std::cout << "===data summary for reduced data with uncertainty:" << std::endl;

    reducedDataSet.AddCellField("entropy", entropyResult);
    reducedDataSet.AddCellField("num_nonzero_prob", numNonZeroProb);
    reducedDataSet.AddCellField("cross_prob", crossProb);

    // reducedDataSet.PrintSummary(std::cout);
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << isovalue;
    std::string isostr = stream.str();

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + "_iso" + isostr + "_" + distribution + "_block" + std::to_string(blocksize) + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(reducedDataSet);

    return 0;
}
