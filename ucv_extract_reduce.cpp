#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>

#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/cont/ArrayHandle.h>

#include <float.h>


int main(int argc, char *argv[])
{
    // init the vtkm (set the backend and log level here)
    vtkm::cont::Initialize(argc, argv);

    if (argc != 3)
    {
        std::cout << "executable <filename> <fieldname>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    std::string fieldName = argv[2];
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

    vtkm::cont::ArrayHandle<vtkm::Id> keyArray;
    vtkm::Id index = 0;

    vtkm::Id blocksize = 4;

    vtkm::Id numberBlockx = xdim / blocksize;
    vtkm::Id numberBlocky = ydim / blocksize;
    vtkm::Id numberBlockz = zdim / blocksize;

    // TODO computing number of coarse grid at each dim
    // computing the id of the coarse grid in the subsequent 
    // three for loop

    for (vtkm::Id k = 0; k < zdim; k++)
    {
        for (vtkm::Id j = 0; j < ydim; j++)
        {
            for (vtkm::Id i = 0; i < xdim; i++)
            {
                //vtkm::Id key = ComputeKey(i, j, k, blocksize);
                // set the key to array
                keyArray.WritePortal.Set(index, key);
                index++;
            }
        }
    }

    return 0;
}
