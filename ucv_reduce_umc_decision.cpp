#include <float.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/DispatcherReduceByKey.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/DispatcherMapTopology.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/ExtractingMinMax.hpp"
#include "ucvworklet/ExtractingMeanStdev.hpp"

#include "ucvworklet/EntropyUniform.hpp"
#include "ucvworklet/EntropyIndependentGaussian.hpp"

int oneDBlocks = 128;
int threadsPerBlock = 128;
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

void initBackend()
{
    // init the vtkh device
    char const *tmp = getenv("UCV_VTKM_BACKEND");
    std::string backend;
    if (tmp == nullptr)
    {
        std::cout << "no UCV_VTKM_BACKEND env, use openmp" << std::endl;
        backend = "openmp";
    }
    else
    {
        backend = std::string(tmp);
    }

    // if (rank == 0)
    //{
    std::cout << "vtkm backend is:" << backend << std::endl;
    //}

    if (backend == "serial")
    {
        vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagSerial());
    }
    else if (backend == "openmp")
    {
        vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP());
    }
    else if (backend == "cuda")
    {
        vtkm::cont::RuntimeDeviceTracker &device_tracker = vtkm::cont::GetRuntimeDeviceTracker();
        device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagCuda());
    }
    else
    {
        std::cerr << " unrecognized backend " << backend << std::endl;
    }
    return;
}

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
    initBackend();
    // init the vtkm (set the backend and log level here)
    vtkm::cont::Initialize(argc, argv);

    if (argc != 5)
    {
        std::cout << "executable <filename> <fieldname> <blocksize> <isovalue>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    int blocksize = std::stoi(argv[3]);
    double isovalue = std::atof(argv[4]);
    std::string distribution = "ig";

    // load the dataset (beetles data set, structured one)
    // TODO, the data set can be distributed between different ranks

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // check the property of the data
    inData.PrintSummary(std::cout);

    // TODO timer start to extract key
    // auto timer1 = std::chrono::steady_clock::now();
    vtkm::cont::Timer timer;
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

    auto coords = inData.GetCoordinateSystem();
    auto bounds = coords.GetBounds();

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
    timer.Stop();
    std::cout << "extractKey time: " << timer.GetElapsedTime() << std::endl;

    timer.Start();

    // indepedent gaussian case
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
    timer.Stop();
    std::cout << "ExtractingMeanStdev time: " << timer.GetElapsedTime() << std::endl;

    timer.Start();
    using WorkletType = EntropyIndependentGaussian;
    using DispatcherEntropyIG = vtkm::worklet::DispatcherMapTopology<WorkletType>;

    DispatcherEntropyIG dispatcherEntropyIG(EntropyIndependentGaussian{isovalue});
    dispatcherEntropyIG.Invoke(reducedDataSet.GetCellSet(), meanArray, stdevArray, crossProb, numNonZeroProb, entropyResult);

    // TODO timer ok to compute uncertianty
    // auto timer4 = std::chrono::steady_clock::now();
    // float EIGaussianTime =
    //    std::chrono::duration<float, std::milli>(timer4 - timer3).count();
    timer.Stop();
    std::cout << "EIGaussianTime time: " << timer.GetElapsedTime() << std::endl;

    // using the same type as the assumption for the output type
    // std::cout << "===data summary for reduced data with uncertainty:" << std::endl;

    reducedDataSet.AddCellField("entropy", entropyResult);
    reducedDataSet.AddCellField("num_nonzero_prob", numNonZeroProb);
    reducedDataSet.AddCellField("cross_prob", crossProb);

    // reducedDataSet.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + "_" + distribution + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(reducedDataSet);



    //using the entropy as the decision resluts
    //to compute another key generation strategy
    //TODO, considering how to map the index between the 
    //reduced data set and the orignal data set


    return 0;
}
