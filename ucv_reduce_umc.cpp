#include <float.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/DispatcherReduceByKey.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

#include <vtkm/worklet/DispatcherMapTopology.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/ExtractingMinMax.hpp"
#include "ucvworklet/ExtractingMeanStdev.hpp"

#include "ucvworklet/EntropyUniform.hpp"
#include "ucvworklet/EntropyIndependentGaussian.hpp"
// compute the entropy based on uniform distribution

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
    vtkm::cont::Initialize(argc, argv);

    if (argc != 6)
    {
        std::cout << "executable <filename> <fieldname> <distribution> <blocksize> <isovalue>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
    std::string distribution = argv[3];
    int blocksize = std::stoi(argv[4]);
    double isovalue = std::atof(argv[5]);

    // load the dataset (beetles data set, structured one)
    // TODO, the data set can be distributed between different ranks

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // check the property of the data
    inData.PrintSummary(std::cout);

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

    // Step3 computing entropy based on reduced data set
    // uniform, indepedent gaussian, multivariant gaussian

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result1;
    vtkm::cont::ArrayHandle<vtkm::Id> result2;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> result3;

    if (distribution == "uni")
    {
        // Step1 creating new key
        vtkm::cont::ArrayHandle<vtkm::Id> keyArrayNew;

        using DispatcherCreateKey = vtkm::worklet::DispatcherMapField<CreateNewKeyWorklet>;
        DispatcherCreateKey dispatcher(CreateNewKeyWorklet{xdim, ydim, zdim,
                                                           numberBlockx, numberBlocky, numberBlockz,
                                                           blocksize});

        dispatcher.Invoke(keyArray, keyArrayNew);


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

        // generate the new data sets with min and max
        // reducedDataSet.AddPointField("ensemble_min", minArray);
        // reducedDataSet.AddPointField("ensemble_max", maxArray);
        // reducedDataSet.PrintSummary(std::cout);

        // output the dataset into the vtk file for results checking
        // std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
        // std::string outputFileName = fileSuffix + std::string("_ReduceDerived.vtk");
        // vtkm::io::VTKDataSetWriter write(outputFileName);
        // write.WriteDataSet(reducedDataSet);

        // uniform distribution
        using WorkletType = EntropyUniform;
        using DispatcherEntropyUniform = vtkm::worklet::DispatcherMapTopology<WorkletType>;

        DispatcherEntropyUniform dispatcherEntropyUniform(EntropyUniform{isovalue});
        dispatcherEntropyUniform.Invoke(reducedDataSet.GetCellSet(), minArray, maxArray, result1, result2, result3);
    }
    else if (distribution == "ig")
    {
        // indepedent gaussian

        // Step1 creating new key
        vtkm::cont::ArrayHandle<vtkm::Id> keyArrayNew;

        using DispatcherCreateKey = vtkm::worklet::DispatcherMapField<CreateNewKeyWorklet>;
        DispatcherCreateKey dispatcher(CreateNewKeyWorklet{xdim, ydim, zdim,
                                                           numberBlockx, numberBlocky, numberBlockz,
                                                           blocksize});

        dispatcher.Invoke(keyArray, keyArrayNew);

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

        using WorkletType = EntropyIndependentGaussian;
        using DispatcherEntropyIG = vtkm::worklet::DispatcherMapTopology<WorkletType>;

        DispatcherEntropyIG dispatcherEntropyIG(EntropyIndependentGaussian{isovalue});
        dispatcherEntropyIG.Invoke(reducedDataSet.GetCellSet(), meanArray, stdevArray, result1, result2, result3);
    }
    else if (distribution == "mg")
    {
        // multivariant gaussian
        // For the multi variant gaussian case, we need to create another key
        // to compute the covariance matrix
        // in order to compute the covariance matrix, the data size
        // assigned to each thread should be 2 times larger than original one
        // assuming each original block is 4*4
        // inorder to compute covariance matrix, we need to access 8*8 data block at once
        // the reduced data is twice smaller as original data
        // since it is based on the large blocks

        // reduced value for uncertainty are per cell
        // for 2d case, each cell contains 4 points 2*2
        // in order to get the value for this 2*2 case
        // the raw data is at least 8*8 and the hixel blcoks are 4

        // two solutions
        // 1 use 8*8 for each thread then compute the u and cov then the final entropy
        // 2 use reduce by key 4*4 to compute u and variance
        // use original data and reduced data to compute the covariance
        // then the cross probability
    }
    else
    {
        throw std::runtime_error("unsupported distribution: " + distribution);
    }

    // using the same type as the assumption for the output type

    std::cout << "===data summary for reduced data with uncertainty:" << std::endl;
    reducedDataSet.AddCellField("cross_prob", result1);
    reducedDataSet.AddCellField("num_nonzero_prob", result2);
    reducedDataSet.AddCellField("entropy", result3);

    reducedDataSet.PrintSummary(std::cout);

    // output the dataset into the vtk file for results checking
    std::string fileSuffix = fileName.substr(0, fileName.size() - 4);
    std::string outputFileName = fileSuffix + "_" + distribution + std::string("_Prob.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.SetFileTypeToBinary();
    write.WriteDataSet(reducedDataSet);

    return 0;
}
