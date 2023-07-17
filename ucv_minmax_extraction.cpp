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
// #include "ucvworklet/MVGaussianWithEnsemble3D.hpp"
#include "ucvworklet/MVGaussianWithEnsemble3DTryLialg.hpp"

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

    if (argc != 5)
    {
        std::cout << "executable [VTK-m options] <filename> <fieldname1> <fieldname2> <blocksize>" << std::endl;
        std::cout << "VTK-m options are:\n";
        std::cout << initResult.Usage << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName1 = argv[2];
    std::string fieldName2 = argv[3];

    int blocksize = std::stoi(argv[4]);

    // load the dataset (beetles data set, structured one)
    // TODO, the data set can be distributed between different ranks

    // create the vtkm data set from the loaded data
    std::cout << "fileName1: " << fileName << std::endl;
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // check the property of the data
    inData.PrintSummary(std::cout);

    // TODO timer start to extract key
    // auto timer1 = std::chrono::steady_clock::now();
    timer.Synchronize();
    timer.Start();

    auto field1 = inData.GetField(fieldName1);
    auto field2 = inData.GetField(fieldName2);

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
    // vtkm::Bounds bounds(0,495,0,831,0,831);
    // supernova
    vtkm::Bounds bounds(-1.02655e+09, 9.42128e+08, -1.02353e+09, 1.00073e+09, -9.16584e+08, 9.32988e+08);

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

    // Step2 extracting ensemble data based on new key
    using DispatcherType = vtkm::worklet::DispatcherReduceByKey<ExtractingMinMax>;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> minArrayF1;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxArrayF1;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> minArrayF2;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> maxArrayF2;
    vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

    auto resolveType = [&](const auto &concrete)
    {
        DispatcherType dispatcher;
        dispatcher.Invoke(keys, concrete, minArrayF1, maxArrayF1);
    };

    field1.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    // generate the new data sets with min and max
    reducedDataSet.AddPointField(fieldName1 + "_ensemble_min", minArrayF1);
    reducedDataSet.AddPointField(fieldName1 + "_ensemble_max", maxArrayF1);

    auto resolveType2 = [&](const auto &concrete)
    {
        DispatcherType dispatcher;
        dispatcher.Invoke(keys, concrete, minArrayF2, maxArrayF2);
    };

    field2.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType2);

    // generate the new data sets with min and max
    reducedDataSet.AddPointField(fieldName2 + "_ensemble_min", minArrayF2);
    reducedDataSet.AddPointField(fieldName2 + "_ensemble_max", maxArrayF2);

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + std::string("_ReduceMinMax_Two_Variables.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(reducedDataSet);

    return 0;
}
