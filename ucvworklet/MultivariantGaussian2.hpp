#ifndef UCV_MULTIVARIANT_GAUSSIAN2_h
#define UCV_MULTIVARIANT_GAUSSIAN2_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>

// this worklet is for the input data that put the different data in a separate array
// for the wind data here https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200
// there are 15 numbers (ensemble extraction) each data is put in a different file

class MultivariantGaussian2 : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MultivariantGaussian2(double isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble, typename OutCellFieldType>

    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutCellFieldType &outCellFieldCProb) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numVertexies = inPointFieldVecEnsemble.GetNumberOfComponents();
        // TODO, extracting data from 4 vertexies and compute the uncertainty things
        // this is supposed to be 4
        //std::cout << "size of numVertexies " << numVertexies << std::endl;
        // this is supposed to be 15
        //std::cout << "ensemble number " << inPointFieldVecEnsemble[0].GetNumberOfComponents() << std::endl;
        

        // this InPointFieldVecEnsemble here is supposed to be the Vec15
        // TODO, compute the mean, cov and cross probability
    
    }

private:
    double m_isovalue;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2_h