#include <vtkm/cont/Initialize.h>
#include <vtkm/io/reader/VTKDataSetReader.h>

int main(int argc, char *argv[])
{
//init the vtkm (set the backend and log level here)
vtkm::cont::Initialize(argc, argv);

if (argc!=2){
    std::cout << "executable <filename>" << std::endl;
    exit(0);
}

std::string fileName = argv[1];
//load the dataset (beetles data set, structured one)
//TODO, the data set can be distributed between different ranks

//create the vtkm data set from the loaded data
vtkm::io::VTKDataSetReader reader(fileName);
vtkm::cont::DataSet inData = reader.ReadDataSet();

//check the property of the data
//the raw_data for testing is 832*832*494
inData.PrintSummary(std::cout);

//TODO try to use the point neighborhood
//skip particular on if it is covered by previous points
//refer to this:
//https://gitlab.kitware.com/vtk/vtk-m/-/blob/release-1.9/vtkm/filter/image_processing/ImageMedian.cxx
//or the convolving small kernels in the vtkm user guide


return 0;
}

