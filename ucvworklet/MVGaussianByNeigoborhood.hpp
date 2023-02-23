#ifndef UCV_MULTIVARIANT_GAUSSIAN_BYNEIGHBOR_h
#define UCV_MULTIVARIANT_GAUSSIAN_BYNEIGHBOR_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/exec/BoundaryState.h>
#include <cmath>
#include "./linalg/ucv_matrix_static_8by8.h"

class MVGaussianByNeigoborhood : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    MVGaussianByNeigoborhood(double isovalue, int numSamples, int blockSize)
        : m_isovalue(isovalue), m_numSamples(numSamples), m_blockSize(blockSize){};

    using ControlSignature = void(CellSetIn,
                                  FieldIn,
                                  FieldOut);

    using ExecutionSignature = void(Boundary, _2, _3);


    template <typename InPointFieldRaw,
              typename OutPointFieldType1>

    VTKM_EXEC void operator()(
        const vtkm::exec::BoundaryState &boundary,
        const InPointFieldRaw &inPointFieldRaw,
        OutPointFieldType1 &outPointMean) const
    {
        // get the ijk
        std::cout << boundary.IJK << std::endl;
    }

private:
    double m_isovalue;
    int m_numSamples;
    int m_blockSize;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN3D_h