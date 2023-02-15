#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>

#include "ContourUncertainEnsemble.h"
#include "ContourUncertainIndependentGaussian.h"
#include "ContourUncertainUniform.h"
#include "SubsampleUncertaintyEnsemble.h"
#include "SubsampleUncertaintyIndependentGaussian.h"
#include "SubsampleUncertaintyUniform.h"

#include <sstream>
#include <iomanip>

int oneDBlocks = 16;
int threadsPerBlock = 16;
#ifdef VTKM_CUDA
// Note: this header will require this file to be compiled with nvcc, but it is required
// for vtkm::cont::cuda::ScheduleParameters.
#include <vtkm/cont/cuda/DeviceAdapterCuda.h>

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
    vtkm::cont::Timer timer{ initResult.Device };
    initBackend(timer);

    std::cout << "initResult.Device: " << initResult.Device.GetName() <<  " timer device: " << timer.GetDevice().GetName() << std::endl;

    if (argc != 6)
    {
        std::cout << "executable [VTK-m options] <filename> <fieldname> <distribution> <blocksize> <isovalue>" << std::endl;
        std::cout << "VTK-m options are:\n";
        std::cout << initResult.Usage << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    std::string distribution = argv[3];
    int blocksize = std::stoi(argv[4]);
    double isovalue = std::atof(argv[5]);

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
    vtkm::cont::DataSet dataset = reader.ReadDataSet();

    // check the property of the data
    dataset.PrintSummary(std::cout);

    // Implementation note: it is typical when running a filter to store the results in
    // a new `DataSet`. However, we are using the same `DataSet` object over and over.
    // This is OK as C++ will properly manage the object and we won't need the data once
    // a filter is run on it. A more important consequence is that when we reuse the
    // `DataSet` object, the old data goes out of scope and any memory it used that is no
    // longer being used gets deleted. This has the desirable side effect of booting
    // data off of a device, which might be important if uniform memory is not being used.

    if (distribution == "uni")
    {
      // uniform
      vtkm::filter::uncertainty::SubsampleUncertaintyUniform subsample;
      subsample.SetBlockSize(blocksize);

      timer.Start();
      dataset = subsample.Execute(dataset);
      timer.Stop();
      std::cout << "extractMinMax time: " << timer.GetElapsedTime() << std::endl;

      vtkm::filter::uncertainty::ContourUncertainUniform contour;
      contour.SetMinField(fieldName + subsample.GetMinSuffix());
      contour.SetMaxField(fieldName + subsample.GetMaxSuffix());
      contour.SetIsoValue(isovalue);

      timer.Start();
      dataset = contour.Execute(dataset);
      timer.Stop();
      std::cout << "EntropyUniformTime time: " << timer.GetElapsedTime() << std::endl;
    }
    else if (distribution == "ig")
    {
      // indepedent gaussian
      vtkm::filter::uncertainty::SubsampleUncertaintyIndependentGaussian subsample;
      subsample.SetBlockSize(blocksize);

      timer.Start();
      dataset = subsample.Execute(dataset);
      timer.Stop();
      std::cout << "ExtractingMeanStdev time: " << timer.GetElapsedTime() << std::endl;

      vtkm::filter::uncertainty::ContourUncertainIndependentGaussian contour;
      contour.SetMeanField(fieldName + subsample.GetMeanSuffix());
      contour.SetStdevField(fieldName + subsample.GetStdevSuffix());
      contour.SetIsoValue(isovalue);

      timer.Start();
      dataset = contour.Execute(dataset);
      timer.Stop();
      std::cout << "EIGaussianTime time: " << timer.GetElapsedTime() << std::endl;
    }
    else if (distribution == "mg")
    {
      // multivariate gaussian
      vtkm::filter::uncertainty::SubsampleUncertaintyEnsemble subsample;
      subsample.SetBlockSize(blocksize);

      timer.Start();
      dataset = subsample.Execute(dataset);
      timer.Stop();
      std::cout << "ExtractingMeanRawTime time: " << timer.GetElapsedTime() << std::endl;

      vtkm::filter::uncertainty::ContourUncertainEnsemble contour;
      contour.SetMeanField(fieldName);
      contour.SetEnsembleField(fieldName + subsample.GetEnsembleSuffix());
      contour.SetIsoValue(isovalue);

      timer.Start();
      dataset = contour.Execute(dataset);
      timer.Stop();
      std::cout << "MVGTime time: " << timer.GetElapsedTime() << std::endl;
    }
    else
    {
        throw std::runtime_error("unsupported distribution: " + distribution);
    }

    // dataset.PrintSummary(std::cout);
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << isovalue;
    std::string isostr = stream.str();

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + "_iso" + isostr + "_" + distribution + "_block" + std::to_string(blocksize) + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(dataset);

    return 0;
}
