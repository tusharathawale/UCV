#include <iostream>
#include <random>
#include <vtkm/Pair.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Algorithm.h>
#include "../worklet/ExtractingMinMax.hpp"

#include "../Fiber.h"
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Timer.h>

int main(int argc, char *argv[])
{
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    std::cout << "initResult.Device: " << initResult.Device.GetName() << std::endl;

    if (argc != 3)
    {
        std::cout << "<executable> <DataFolder> <NumEns>" << std::endl;
        exit(0);
    }

    std::string dataFolder = std::string(argv[1]);
    int NumEns = std::stoi(argv[2]);

    // compute the min and max through the mean+-stdev for two variables
    std::string CurlField = "curlZ";
    std::string VorField = "vorticityMagnitude";

    // the name of the file is mistyped, the value is actually the curl
    std::string MeanCurlFile = dataFolder + "/curlZ/meanVol/meanVorticity.vtk";
    std::string DevDevFile = dataFolder + "/curlZ/devVol/devVorticity.vtk";

    std::string MeanVorFile = dataFolder + "/vorticityMagnitude/meanVol/meanVorticity.vtk";
    std::string DevVorFile = dataFolder + "/vorticityMagnitude/devVol/devVorticity.vtk";

    // get mean for the curl
    vtkm::io::VTKDataSetReader MeanCurlReader(MeanCurlFile);
    vtkm::cont::DataSet MeanCurlData = MeanCurlReader.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> MeanCurlDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(inData.GetField("meanVorticity").GetData(), MeanCurlDataArray);

    // get the cellset
    auto cellSet = MeanCurlData.GetCellSet();
    vtkm::cont::CellSetStructured<3> structCellSet =
        cellSet.AsCellSet<vtkm::cont::CellSetStructured<3>>();

    vtkm::Id3 pointDims = structCellSet.GetPointDimensions();
    std::cout << "------" << std::endl;
    std::cout << "point dim: " << pointDims[0] << " " << pointDims[1] << " " << pointDims[2] << std::endl;

    // get dev for the curl
    vtkm::io::VTKDataSetReader DevCurlReader(MeanCurlFile);
    vtkm::cont::DataSet DevCurlData = DevCurlReader.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> MeanDevDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(DevCurlData.GetField("devVorticity").GetData(), MeanDevDataArray);

    // get mean for the vorticity

    // get dev for the vorticity
}