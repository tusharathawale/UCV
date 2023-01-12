
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/worklet/DispatcherReduceByKey.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/Initialize.h>

#include "ucvworklet/CreateNewKey.hpp"
#include "ucvworklet/MVGaussianWithEnsemble2D.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using SupportedTypesVec = vtkm::List<vtkm::Vec<double, 15>>;

int main(int argc, char *argv[])
{
    vtkm::cont::Initialize(argc, argv);
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

    std::vector<vtkm::cont::ArrayHandle<vtkm::Float64>> componentArrays;

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

        // ArrayHandleSOA will group the components in each array and
        // store them into Vec15 automatically.
        // transfer vector into tha array handle
        componentArrays.push_back(vtkm::cont::make_ArrayHandleMove(std::move(d)));
    }

    using Vec15 = vtkm::Vec<double, 15>;
    vtkm::cont::ArrayHandleSOA<Vec15> soaArray(componentArrays);

    vtkmDataSet.AddPointField("ensembles", soaArray);

    // check results

    // std::string outputFileName = "./wind_pressure_200.vtk";
    // vtkm::io::VTKDataSetWriter write(outputFileName);
    // write.WriteDataSet(vtkmDataSet);

    // let the data set go through the multivariant gaussian filter
    using WorkletType = MVGaussianWithEnsemble2D;
    using DispatcherType = vtkm::worklet::DispatcherMapTopology<WorkletType>;

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> crossProbability;
    double isovalue = 0.3;
    auto resolveType = [&](const auto &concrete)
    {
        DispatcherType dispatcher(MVGaussianWithEnsemble2D{isovalue});
        dispatcher.Invoke(vtkmDataSet.GetCellSet(), concrete, crossProbability);
    };

    vtkmDataSet.GetField("ensembles").GetData().CastAndCallForTypes<SupportedTypesVec, VTKM_DEFAULT_STORAGE_LIST>(resolveType);

    // check results
    vtkmDataSet.AddCellField("cross_prob", crossProbability);
    vtkmDataSet.PrintSummary(std::cout);
    std::string outputFileName = "./wind_pressure_200.vtk";
    vtkm::io::VTKDataSetWriter write(outputFileName);
    write.WriteDataSet(vtkmDataSet);
}