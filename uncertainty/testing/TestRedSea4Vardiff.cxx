#include <iostream>
#include <random>
#include <vtkm/Pair.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Algorithm.h>
#include "../worklet/ExtractingMinMaxFromMeanDev.hpp"
#include "../worklet/ComputeDiff.hpp"


#include "../Fiber4Var.h"
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Timer.h>

int main(int argc, char *argv[])
{
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    std::cout << "initResult.Device: " << initResult.Device.GetName() << std::endl;

    if (argc != 4)
    {
        std::cout << "<executable> <DataFolder> <Approach> <NumSamples>" << std::endl;
        exit(0);
    }

    std::string dataFolder = std::string(argv[1]);
    //int NumEns = 20;

    std::string Approach = std::string(argv[2]);
    int NumSamples = std::stoi(argv[3]); // this only work when appraoch is MonteCarlo

    if (Approach == "MonteCarlo" || Approach == "ClosedForm" || Approach == "Mean")
    {
    }
    else
    {
        std::cout << "Approach should be MonteCarlo or ClosedFrom" << std::endl;
        exit(0);
    }

    // compute the min and max through the mean+-stdev for two variables
    
    std::cout << "got here 1" << std::endl;
    // get mean for the curl
    vtkm::io::VTKDataSetReader ClosedForm(dataFolder + "redSea3VarOutputClosedForm666.vtk");
    vtkm::cont::DataSet ClosedFormData = ClosedForm.ReadDataSet();

    vtkm::cont::ArrayHandle<vtkm::FloatDefault> ClosedFormArray;
    vtkm::cont::ArrayCopyShallowIfPossible(ClosedFormData.GetField("ClosedForm").GetData(), ClosedFormArray);

    // get the cellset

    std::cout << "got here 2" << std::endl;
    // get dev for the curl
    vtkm::io::VTKDataSetReader MonteCarlo(dataFolder + "redSea3VarOutputMonteCarlo1000.vtk");
    vtkm::cont::DataSet MonteCarloData = MonteCarlo.ReadDataSet();
 
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> MonteCarloArray;
    vtkm::cont::ArrayCopyShallowIfPossible(MonteCarloData.GetField("MonteCarlo").GetData(), MonteCarloArray);


    vtkm::cont::Invoker invoke;

    //array for diff feild
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> diffFeild;

    invoke(ComputeDiffSquare{}, ClosedFormArray, MonteCarloArray, diffFeild);

    
    //find the max value

    vtkm::FloatDefault max = 0;

    vtkm::FloatDefault sum = 0.0;
    //std::cout << "diffFeild.GetNumberOfValues(): " << diffFeild.GetNumberOfValues() << std::endl;
    for (vtkm::Id i = 0; i < diffFeild.GetNumberOfValues(); i++)
    {
        
        if(diffFeild.ReadPortal().Get(i) > max)
        {
            max = diffFeild.ReadPortal().Get(i);
        }
        
        sum += diffFeild.ReadPortal().Get(i);
        //std::cout << sum << std::endl;
    }
    vtkm::FloatDefault total = sum / diffFeild.GetNumberOfValues();
    std::cout << "totalsum: " << total << std::endl;
    std::cout << "max: " << max << std::endl;

    
}