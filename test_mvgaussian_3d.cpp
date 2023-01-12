
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/Initialize.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/MVGaussianWithEnsemble3D.hpp"
#include "ucvworklet/ExtractingMeanRaw.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32,
                                  vtkm::Id>;

void testMVGaussian3d()
{
    // create data sets with 8*8 dims for testing
    vtkm::Id xdim = 8;
    vtkm::Id ydim = 8;
    vtkm::Id zdim = 8;

    const vtkm::Id3 dims(xdim, ydim, zdim);
    vtkm::cont::DataSetBuilderUniform dataSetBuilder;

    vtkm::cont::DataSet testDataSet = dataSetBuilder.Create(dims);

    testDataSet.PrintSummary(std::cout);

    // use multivariant gaussian for this dataset

    const int blocksize = 4;

    if (xdim % blocksize != 0 || ydim % blocksize != 0 || zdim % blocksize != 0)
    {
        // if the data size is not divided by blocksize
        // we can reample or padding the data set before hand
        throw std::runtime_error("only support the case where xyz dim is diveide dy blocksize");
    }

    vtkm::Id numberBlockx = xdim % blocksize == 0 ? xdim / blocksize : xdim / blocksize + 1;
    vtkm::Id numberBlocky = ydim % blocksize == 0 ? ydim / blocksize : ydim / blocksize + 1;
    vtkm::Id numberBlockz = zdim % blocksize == 0 ? zdim / blocksize : zdim / blocksize + 1;

    auto keyArray =
        vtkm::cont::ArrayHandleCounting<vtkm::Id>(0, 1, static_cast<vtkm::Id>(xdim * ydim * zdim));

    testDataSet.AddPointField("data", keyArray);

    vtkm::cont::ArrayHandle<vtkm::Id> keyArrayNew;
    using DispatcherCreateKey = vtkm::worklet::DispatcherMapField<CreateNewKeyWorklet>;
    DispatcherCreateKey dispatcher(CreateNewKeyWorklet{xdim, ydim, zdim,
                                                       numberBlockx, numberBlocky, numberBlockz,
                                                       blocksize});

    dispatcher.Invoke(keyArray, keyArrayNew);

    testDataSet.AddPointField("keyarray", keyArrayNew);
    // std::cout <<"summary new array" << std::endl;
    // printSummary_ArrayHandle(keyArrayNew, std::cout);

    std::string outputFileName = std::string("TestMGData.vtk");
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(testDataSet);

    // send data into the worklet to compute the cross probabilities
    using WorkletType = ExtractingMeanRaw;
    using Dispatcher = vtkm::worklet::DispatcherReduceByKey<WorkletType>;

    // the value here should be same with the elements in each hixel
    using VecType = vtkm::Vec<vtkm::FloatDefault, blocksize * blocksize * blocksize>;
    vtkm::cont::ArrayHandle<VecType> SOARawArray;
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> meanArray;
    // Pay attention to transfer the arrayHandle into the Keys type
    vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

    auto resolveType = [&](const auto &concrete)
    {
        Dispatcher dispatcher;
        dispatcher.Invoke(keys, concrete, meanArray, SOARawArray);
    };

    testDataSet.GetField("data").GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    std::cout << "SOARawArray summary:" << std::endl;
    vtkm::cont::printSummary_ArrayHandle(SOARawArray, std::cout);
    std::cout << "meanArray summary:" << std::endl;
    vtkm::cont::printSummary_ArrayHandle(meanArray, std::cout);

    // Then using the MVGaussianWithEnsemble3D to extract the cross probability
    // create reduced dataset
    // using the ensemble filter for the reduced data set
    const vtkm::Id3 reducedDims(numberBlockx, numberBlocky, numberBlockz);
    vtkm::cont::DataSet reducedDataSet = dataSetBuilder.Create(reducedDims);


    using WorkletTypeMVG = MVGaussianWithEnsemble3D;
    using DispatcherTypeMVG = vtkm::worklet::DispatcherMapTopology<WorkletTypeMVG>;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProbability;
    double isovalue = 100;

    DispatcherTypeMVG dispatcherMVG(MVGaussianWithEnsemble3D{isovalue, 1000});
    dispatcherMVG.Invoke(reducedDataSet.GetCellSet(), SOARawArray, meanArray, crossProbability);

    std::cout << "crossProbability summary:" << std::endl;
    vtkm::cont::printSummary_ArrayHandle(crossProbability, std::cout);
    
}

int main(int argc, char *argv[])
{
    vtkm::cont::Initialize(argc, argv);
    testMVGaussian3d();
}