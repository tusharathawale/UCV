#ifndef UCV_EXTRACTING_BYNEIGHBOR_ENS_h
#define UCV_EXTRACTING_BYNEIGHBOR_ENS_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/exec/BoundaryState.h>

// extracting specific ensemble element
class ExtractingByNeigoborhoodEnsTwoFields : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    ExtractingByNeigoborhoodEnsTwoFields(vtkm::IdComponent blockSize,
                                vtkm::Id xdim, vtkm::Id ydim, vtkm::Id zdim, vtkm::Id ensId)
        : m_blockSize(blockSize),
          m_xdim(xdim), m_ydim(ydim), m_zdim(zdim), m_ensId(ensId){};

    // the second f
    using ControlSignature = void(CellSetIn,
                                  WholeArrayIn,
                                  FieldOut,
                                  WholeArrayIn,
                                  FieldOut);
    // workindex will decrease the execution speed of function call
    using ExecutionSignature = void(Boundary, _2, _3, _4, _5);

    template <typename InPointFieldRaw,
              typename OutputType>
    VTKM_EXEC void operator()(
        const vtkm::exec::BoundaryState &boundary,
        const InPointFieldRaw &inPointFieldPortal1,
        OutputType &ensValue1,
        const InPointFieldRaw &inPointFieldPortal2,
        OutputType &ensValue2) const
    {
        // get the ijk (the simple map field can not get the ijk)
        // std::cout << boundary.IJK << std::endl;

        // get the original bounding box
        // left and right for bounding box for each dim
        vtkm::Id l_i = boundary.IJK[0] * m_blockSize;
        vtkm::Id l_j = boundary.IJK[1] * m_blockSize;
        vtkm::Id l_k = boundary.IJK[2] * m_blockSize;

        vtkm::Id r_i = (boundary.IJK[0] + 1) * m_blockSize;
        vtkm::Id r_j = (boundary.IJK[1] + 1) * m_blockSize;
        vtkm::Id r_k = (boundary.IJK[2] + 1) * m_blockSize;

        vtkm::Id ensIndex = 0;
        for (vtkm::Id k = l_k; k < r_k; k++)
        {
            for (vtkm::Id j = l_j; j < r_j; j++)
            {
                for (vtkm::Id i = l_i; i < r_i; i++)
                {
                    vtkm::Id index = k * m_xdim * m_ydim + j * m_xdim + i;
                    if (ensIndex == this->m_ensId)
                    {
                        ensValue1 = inPointFieldPortal1.Get(index);
                        ensValue2 = inPointFieldPortal2.Get(index);
                        return;
                    }
                    ensIndex++;
                }
            }
        }
    }

private:
    vtkm::IdComponent m_blockSize;
    vtkm::Id m_xdim;
    vtkm::Id m_ydim;
    vtkm::Id m_zdim;
    vtkm::Id m_ensId;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN3D_h
