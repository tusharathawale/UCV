
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/ArrayHandleSOA.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/MultivariantGaussian.hpp"
#include "ucvworklet/MultivariantGaussian2.hpp"

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
void testMVGaussian1()
{
    // create data sets with 8*8 dims for testing
    vtkm::Id xdim = 9;
    vtkm::Id ydim = 9;
    vtkm::Id zdim = 1;

    const vtkm::Id3 dims(xdim, ydim, zdim);
    vtkm::cont::DataSetBuilderUniform dataSetBuilder;

    vtkm::cont::DataSet testDataSet = dataSetBuilder.Create(dims);

    testDataSet.PrintSummary(std::cout);

    // use multivariant gaussian for this dataset

    int blocksize = 8;

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
    using WorkletType = MultivariantGaussian;
    using Dispatcher = vtkm::worklet::DispatcherReduceByKey<WorkletType>;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProbability;
    vtkm::worklet::Keys<vtkm::Id> keys(keyArrayNew);

    auto resolveType = [&](const auto &concrete)
    {
        Dispatcher dispatcher(MultivariantGaussian{1, blocksize});
        dispatcher.Invoke(keys, concrete, crossProbability);
    };

    testDataSet.GetField("data").GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    std::cout << "crossProbability:" <<std::endl;
    printSummary_ArrayHandle(crossProbability, std::cout);

}

// dir should contains the "/"
void testMVGaussian2()
{
    // assuming the ensemble data set is already been extracted out
    // we test results by this dataset
    // https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200
    // to make sure the reuslts are correct

    // load data set, the dim is m*n and for each point there are k ensemble values
    // the data set comes from here https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200

    vtkm::Id xdim = 240;
    vtkm::Id ydim = 121;
    vtkm::Id zdim = 1;

    const vtkm::Id3 dims(xdim, ydim, zdim);
    vtkm::cont::DataSetBuilderUniform dataSetBuilder;

    vtkm::cont::DataSet vtkmDataSet = dataSetBuilder.Create(dims);

    std::string windDataDir = "./wind_pressure_200/";

    using Vec15 = vtkm::Vec<double, 15>;
    vtkm::cont::ArrayHandleSOA<Vec15> soaArray;

    for (vtkm::IdComponent i = 0; i < 15; i++)
    {
        std::string filename = windDataDir + "Lead_33_" + std::to_string(i) + ".txt";

        // load the data
        std::ifstream fin(filename.c_str());

        std::vector<double> d;
        float element;
        while (fin >> element)
        {
            d.push_back(element);
        }

        // TODO, how to set the data to the soaArray?

    }

    vtkmDataSet.AddPointField("ensembles", soaArray);
    vtkmDataSet.PrintSummary(std::cout);

    // check results

    // std::string outputFileName = "./wind_pressure_200.vtk";
    // vtkm::io::VTKDataSetWriter write(outputFileName);
    // write.WriteDataSet(vtkmDataSet);

    // let the data set go through the multivariant gaussian filter
    using WorkletType = MultivariantGaussian2;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProbability;

    // DispatcherType dispatcher(MultivariantGaussian2{isovalue});
    // dispatcher.Invoke(vtkmDataSet.GetCellSet(), minArray, maxArray, result1, result2, result3);
}

int main(int argc, char *argv[])
{

    testMVGaussian1();
    //testMVGaussian2();
}