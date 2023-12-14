#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include "vtkmlib/DataSetConverters.h"

int main(int argc, char *argv[])
{

    // load vtk data
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0]
                  << " Filename" << std::endl;
        return EXIT_FAILURE;
    }

    std::string filename = argv[1];

    // Read all the data from the file
    vtkSmartPointer<vtkXMLImageDataReader> reader =
        vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // get the specific unstructureGridData and check the results
    auto readDataSet = reader->GetOutput();

    vtkDataSet *input = vtkDataSet::SafeDownCast(readDataSet);
    input->Print(std::cout);

    // transfer it to vktm data
    vtkm::cont::DataSet in = tovtkm::Convert(input, tovtkm::FieldsFlag::PointsAndCells);

    std::cout << "checking input vtkm data" << std::endl;
    in.PrintSummary(std::cout);
}