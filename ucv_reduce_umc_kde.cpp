#include <float.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapField.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/Timer.h>



#include "ucvworklet/CreateNewKey.hpp"


#include "ucvworklet/HelperProbKDE.hpp"
#include "ucvworklet/EntropyKDE.hpp"
#include "ucvworklet/ExtractingRaw.hpp"

#include <sstream>
#include <iomanip>

//
// #endif // VTKM_CUDA


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
        std::cout << "executable [VTK-m options] <filename> <fieldname> <distribution> <isovalue>" << std::endl;
        std::cout << "VTK-m options are:\n";
        std::cout << initResult.Usage << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    std::string distribution = argv[3];
    double isovalue = std::atof(argv[4]);
    const int blocksize  = 4 ;

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
    //vtkm::Bounds bounds(-1.02655e+09, 9.42128e+08, -1.02353e+09, 1.00073e+09, -9.16584e+08, 9.32988e+08);

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

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> postiveProb;


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

    if (distribution == "kde")
    {

        using WorkletType = ExtractingRaw;
        using DispatcherType = vtkm::worklet::DispatcherReduceByKey<WorkletType>;

        using VecType = vtkm::Vec<vtkm::FloatDefault, blocksize * blocksize * blocksize>;
        vtkm::cont::ArrayHandle<VecType> SOARawArray;
        // Pay attention to transfer the arrayHandle into the Keys type
        vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

        auto resolveType = [&](const auto &concrete)
        {
            DispatcherType dispatcher;
            dispatcher.Invoke(keys, concrete, SOARawArray);
        };

        field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
            resolveType);

        // TODO timer ok to extract property
        // auto timer3 = std::chrono::steady_clock::now();
        // float extractMinMax =
        //    std::chrono::duration<float, std::milli>(timer3 - timer2).count();
        timer.Synchronize();
        timer.Stop();
        std::cout << "ExtractingRaw time: " << timer.GetElapsedTime() * 1000 << std::endl;
        timer.Synchronize();
        timer.Start();
        // generate the new data sets with min and max
        reducedDataSet.AddPointField("ensemble_raw", SOARawArray);
        // reducedDataSet.PrintSummary(std::cout);

        // output the dataset into the vtk file for results checking
        //std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
        //std::string outputFileName = fileSuffix + "_" + fieldName + std::string("_ReduceMinMax.vtk");
        //vtkm::io::VTKDataSetWriter write(outputFileName);
        //write.WriteDataSet(reducedDataSet);

        // compute kde distribution
        using DispatcherProbKDE = vtkm::worklet::DispatcherMapField<HelperProbKDE>;

        DispatcherProbKDE dispatcherProbKDE(HelperProbKDE{isovalue});
        dispatcherProbKDE.Invoke(SOARawArray,postiveProb);
        
        // compute entropy things
        // go through cell by points
        using DispatcherEntropyKDE = vtkm::worklet::DispatcherMapTopology<EntropyKDE>;

        DispatcherEntropyKDE dispatcherEntropyKDE(EntropyKDE{isovalue});
        dispatcherEntropyKDE.Invoke(reducedDataSet.GetCellSet(), postiveProb, crossProb, numNonZeroProb, entropyResult);

        //reducedDataSet.AddPointField("postiveProb", postiveProb);
        //reducedDataSet.PrintSummary(std::cout);

        timer.Synchronize();
        timer.Stop();
        std::cout << "EntropyKDE time: " << timer.GetElapsedTime() * 1000 << std::endl;
    }
    else
    {
        std::cout << "We only support KDE for this case" << std::endl;
        exit(0);
    }

    // using the same type as the assumption for the output type
    // std::cout << "===data summary for reduced data with uncertainty:" << std::endl;

    reducedDataSet.AddCellField("entropy", entropyResult);
    reducedDataSet.AddCellField("num_nonzero_prob", numNonZeroProb);
    reducedDataSet.AddCellField("cross_prob", crossProb);

    reducedDataSet.PrintSummary(std::cout);
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << isovalue;
    std::string isostr = stream.str();

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + "_kde_iso" + isostr + "_" + distribution + "_block" + std::to_string(blocksize) + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(reducedDataSet);

    return 0;
}
