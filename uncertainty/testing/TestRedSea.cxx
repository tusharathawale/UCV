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

    if (argc != 2)
    {
        std::cout << "<executable> <DataFolder>" << std::endl;
        exit(0);
    }

    std::string dataFolder = std::string(argv[1]);
    int NumEns = 20;

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
    vtkm::cont::ArrayCopyShallowIfPossible(MeanCurlData.GetField("meanVorticity").GetData(), MeanCurlDataArray);

    // get the cellset
    auto cellSet = MeanCurlData.GetCellSet();
    vtkm::cont::CellSetStructured<3> structCellSet =
        cellSet.AsCellSet<vtkm::cont::CellSetStructured<3>>();

    vtkm::Id3 pointDims = structCellSet.GetPointDimensions();
    std::cout << "------" << std::endl;
    std::cout << "point dim: " << pointDims[0] << " " << pointDims[1] << " " << pointDims[2] << std::endl;

    // get dev for the curl
    vtkm::io::VTKDataSetReader DevCurlReader(DevDevFile);
    vtkm::cont::DataSet DevCurlData = DevCurlReader.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> DevCurlDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(DevCurlData.GetField("devVorticity").GetData(), DevCurlDataArray);

    // get mean for the vorticity
    vtkm::io::VTKDataSetReader MeanVorReader(MeanVorFile);
    vtkm::cont::DataSet MeanVorData = MeanVorReader.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> MeanVorDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(MeanVorData.GetField("meanVorticity").GetData(), MeanVorDataArray);

    // get dev for the vorticity
    vtkm::io::VTKDataSetReader DevVorReader(DevVorFile);
    vtkm::cont::DataSet DevVorData = MeanVorReader.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> DevVorDataArray;
    vtkm::cont::ArrayCopyShallowIfPossible(DevVorData.GetField("devVorticity").GetData(), DevVorDataArray);
}