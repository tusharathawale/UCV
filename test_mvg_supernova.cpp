
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/Initialize.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/MVGaussianWithEnsemble3DTryLialg.hpp"
#include "ucvworklet/ExtractingMeanRaw.hpp"

#include <vtkm/cont/Timer.h>
#include <mpi.h>

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32,
                                  vtkm::Id>;

std::string backend = "openmp";

void initBackend(vtkm::cont::Timer &timer, int rank)
{
    // init the vtkh device
    char const *tmp = getenv("UCV_VTKM_BACKEND");

    if (tmp == nullptr)
    {
        std::cout << "no UCV_VTKM_BACKEND env, use openmp" << std::endl;
        backend = "openmp";
    }
    else
    {
        backend = std::string(tmp);
    }

    if (rank == 0)
    {
        std::cout << "vtkm backend is:" << backend << std::endl;
    }

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

void callWorklet(vtkm::cont::DataSet inData, double iso, int numSamples, int blocksize, int rank, int blockid, std::string outputStr)
{
    // 3d structured case for multi variant gaussian
    std::string fieldName = "Iron";

    // Step 0 declaring necessary field and data set
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

    // go through all points and set the specific key
    vtkm::Id xdim = pointDims[0];
    vtkm::Id ydim = pointDims[1];
    vtkm::Id zdim = pointDims[2];

    auto keyArray =
        vtkm::cont::ArrayHandleCounting<vtkm::Id>(0, 1, static_cast<vtkm::Id>(xdim * ydim * zdim));

    vtkm::Id numberBlockx = xdim % blocksize == 0 ? xdim / blocksize : xdim / blocksize + 1;
    vtkm::Id numberBlocky = ydim % blocksize == 0 ? ydim / blocksize : ydim / blocksize + 1;
    vtkm::Id numberBlockz = zdim % blocksize == 0 ? zdim / blocksize : zdim / blocksize + 1;

    const vtkm::Id3 reducedDims(numberBlockx, numberBlocky, numberBlockz);

    auto coords = inData.GetCoordinateSystem();
    auto bounds = coords.GetBounds();

    auto reducedOrigin = bounds.MinCorner();

    vtkm::FloatDefault spacex = (bounds.X.Max - bounds.X.Min) / (numberBlockx - 1);
    vtkm::FloatDefault spacey = (bounds.Y.Max - bounds.Y.Min) / (numberBlocky - 1);
    vtkm::FloatDefault spacez = (bounds.Z.Max - bounds.Z.Min) / (numberBlockz - 1);

    vtkm::Vec3f_64 reducedSpaceing(spacex, spacey, spacez);

    vtkm::cont::DataSetBuilderUniform dataSetBuilder;

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

    // Step2
    //   extracting the mean and rawdata for each hixel block
    //   the raw data is used to compute the covariance matrix
    if (xdim % 4 != 0 || ydim % 4 != 0 || zdim % 4 != 0)
    {
        // if the data size is not divided by blocksize
        // we can reample or padding the data set before hand
        // it will be convenient to compute cov matrix by this way
        throw std::runtime_error("only support blocksize = 4 and the case where xyz dim is diveide dy blocksize for current mg");
    }

    // the value here should be same with the elements in each hixel
    using WorkletType = ExtractingMeanRaw;
    using DispatcherType = vtkm::worklet::DispatcherReduceByKey<WorkletType>;

    // this should be modified if the blocksize change
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

    // Step3 computing the cross probability
    using WorkletTypeMVG = MVGaussianWithEnsemble3DTryLialg;
    using DispatcherTypeMVG = vtkm::worklet::DispatcherMapTopology<WorkletTypeMVG>;

    DispatcherTypeMVG dispatcherMVG(WorkletTypeMVG{iso, numSamples});
    dispatcherMVG.Invoke(reducedDataSet.GetCellSet(), SOARawArray, meanArray, crossProb, numNonZeroProb, entropyResult);

    // https://public.kitware.com/pipermail/paraview/2005-November/002188.html
    // output the data for testing
    // multiple files
    // create the vtm file at last
    if (outputStr == "true")
    {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << iso;
        std::string isostr = stream.str();

        reducedDataSet.AddCellField("cross_prob_" + isostr, crossProb);
        reducedDataSet.AddCellField("num_nonzero_prob_" + isostr, numNonZeroProb);
        reducedDataSet.AddCellField("entropy_" + isostr, entropyResult);

        std::string outputFileName = "./supernova_reduced/output_" + std::to_string(blockid) + ".vtk";
        vtkm::io::VTKDataSetWriter writeCross(outputFileName);
        writeCross.WriteDataSet(reducedDataSet);
        std::cout << "ok to write data to " << outputFileName << std::endl;
    }
}

int main(int argc, char *argv[])
{

    int rc = MPI_Init(&argc, &argv);
    int rank;
    int numProcesses;
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // detect vtkm related commnad line parameters
    vtkm::cont::Initialize(argc, argv);

    if (argc != 5)
    {
        if (rank == 0)
        {
            std::cout << "<executable> <iso> <blocksize> <num of samples> <output>" << std::endl;
        }
        exit(0);
    }
    // 0.3 for the iron value can be a good option
    double isovalue = std::stod(argv[1]);
    int blocksize = std::stoi(argv[2]);
    int num_samples = std::stoi(argv[3]);
    std::string outputStr = argv[4];

    if (rank == 0)
    {
        std::cout << "iso value is: " << isovalue << " blocksize: " << blocksize << " num_samples is: " << num_samples << std::endl;
    }

    vtkm::cont::Timer timer;
    initBackend(timer, rank);
    if (rank == 0)
    {
        std::cout << "timer device: " << timer.GetDevice().GetName() << std::endl;
    }

    int numDataBlocks = 512;
    // this is for testing in small scale
    // int numDataBlocks = 128;

    if (numProcesses > numDataBlocks)
    {
        if (rank == 0)
        {
            std::cout << "number of processes should <=512" << std::endl;
        }
    }

    std::unordered_map<int, vtkm::cont::DataSet> dsMap;
    // load data
    std::string dataDir = "./supernova_decompose";
    for (int blockid = 0; blockid < numDataBlocks; blockid++)
    {
        if (blockid % numProcesses == rank)
        {
            // current rank load the ith slice
            char str[32];
            sprintf(str, "%03d", blockid);
            std::string fileName = dataDir + "/output." + std::string(str) + ".vtk";
            std::cout << "rank " << rank << " load file " << fileName << std::endl;
            vtkm::io::VTKDataSetReader reader(fileName);
            vtkm::cont::DataSet inData = reader.ReadDataSet();
            dsMap[blockid] = inData;
        }
    }

    // process data for each rank
    // do the processing for each member in the list
    // time it
    MPI_Barrier(MPI_COMM_WORLD);
    timer.Start();

    for (auto it = dsMap.begin(); it != dsMap.end(); it++)
    {
        callWorklet(it->second, isovalue, num_samples, blocksize, rank, it->first, outputStr);
    }

    // maybe add more operations here
    // such as adding the collective operations based on entropy
    // and compute some results back
    // do some further analysis here

    MPI_Barrier(MPI_COMM_WORLD);
    timer.Stop();

    if (rank == 0 && outputStr == "true")
    {
        // write the metadata file for visit cases
        std::ofstream metafile;
        metafile.open("supernova_reduced.visit");
        metafile << "!NBLOCKS " << numDataBlocks << std::endl;
        for (int i = 0; i < numDataBlocks; i++)
        {
            std::string outputFileName = "./supernova_reduced/output_" + std::to_string(i) + ".vtk\n";
            metafile << outputFileName;
        }
        metafile.close();
    }

    if (rank == 0)
    {
        std::cout << "execution time for rank 0: " << timer.GetElapsedTime() * 1000 << std::endl;
    }

    MPI_Finalize();

    return 0;
}
