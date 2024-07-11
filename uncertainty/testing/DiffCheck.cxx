#include <iostream>
#include <random>
#include <vtkm/Pair.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Algorithm.h>
#include "../worklet/ExtractingMinMaxFromMeanDev.hpp"
#include "../worklet/ComputeDiff.hpp"

#include "../Fiber3Var.h"
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/Timer.h>

int main(int argc, char *argv[])
{
    vtkm::cont::InitializeResult initResult = vtkm::cont::Initialize(
        argc, argv, vtkm::cont::InitializeOptions::DefaultAnyDevice);
    std::cout << "initResult.Device: " << initResult.Device.GetName() << std::endl;

    if (argc != 6)
    {
        std::cout << "<executable> <DataFolder> <File1> <Feild1> <File2> <Field2>" << std::endl;
        exit(0);
    }

    std::string dataFolder = std::string(argv[1]);
    //int NumEns = 20;

    std::string File1 = std::string(argv[2]);
    std::string Field1 = std::string(argv[3]);
    std::string File2 = std::string(argv[4]);
    std::string Field2 = std::string(argv[5]);


    // compute the min and max through the mean+-stdev for two variables
    
    
    // get mean for the curl
    vtkm::io::VTKDataSetReader ClosedForm(dataFolder + File1);
    vtkm::cont::DataSet ClosedFormData = ClosedForm.ReadDataSet();
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> ClosedFormArray;
    vtkm::cont::ArrayCopyShallowIfPossible(ClosedFormData.GetField(Field1).GetData(), ClosedFormArray);
    std::cout << "Read " + Field1 + " From " + File1 << std::endl;

    // get the cellset

    
    // get dev for the curl
    vtkm::io::VTKDataSetReader MonteCarlo(dataFolder + File2);
    vtkm::cont::DataSet MonteCarloData = MonteCarlo.ReadDataSet();
 
    vtkm::cont::ArrayHandle<vtkm::FloatDefault> MonteCarloArray;
    vtkm::cont::ArrayCopyShallowIfPossible(MonteCarloData.GetField(Field2).GetData(), MonteCarloArray);
    std::cout << "Read " + Field2 + " From " + File2 << std::endl;

    vtkm::cont::Invoker invoke;

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
        //checks for any nan values in either of the fields
        if (std::isnan(ClosedFormArray.ReadPortal().Get(i)))
        {
            std::cout << "nan in File1 value found at: " << i << std::endl;
        }
        else if (std::isnan(MonteCarloArray.ReadPortal().Get(i)))
        {
            std::cout << "nan in File2  found at: " << i << std::endl;
        }else if (std::isnan(diffFeild.ReadPortal().Get(i)))
        {
            std::cout << "nan in Diff  found at: " << i << std::endl;
        }
        
        
    }
    vtkm::FloatDefault total = sum / diffFeild.GetNumberOfValues();
    std::cout << "total sum: " << total << std::endl;
    std::cout << "max: " << max << std::endl;

    
}