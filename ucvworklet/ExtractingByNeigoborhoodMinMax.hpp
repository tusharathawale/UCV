#ifndef UCV_EXTRACTING_BYNEIGHBOR_MIN_MAX_h
#define UCV_EXTRACTING_BYNEIGHBOR_MIN_MAX_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/exec/BoundaryState.h>

class ExtractingByNeigoborhoodMinMax : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    ExtractingByNeigoborhoodMinMax(vtkm::IdComponent blockSize,
                                   vtkm::Id xdim, vtkm::Id ydim, vtkm::Id zdim)
        : m_blockSize(blockSize),
          m_xdim(xdim), m_ydim(ydim), m_zdim(zdim){};

    // the second f
    using ControlSignature = void(CellSetIn,
                                  WholeArrayIn,
                                  FieldOut,
                                  FieldOut);
    //workindex will decrease the execution speed of function call
    using ExecutionSignature = void(Boundary, _2, _3, _4);

    template <typename InPointFieldRaw,
              typename OutputType>

    VTKM_EXEC void operator()(
        const vtkm::exec::BoundaryState &boundary,
        const InPointFieldRaw &inPointFieldPortal,
        OutputType &minValue, OutputType &maxValue) const
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

        OutputType boxMin = std::numeric_limits<OutputType>::max();
        vtkm::FloatDefault boxMax = std::numeric_limits<OutputType>::min();

        for (vtkm::Id k = l_k; k < r_k; k++)
        {
            for (vtkm::Id j = l_j; j < r_j; j++)
            {
                for (vtkm::Id i = l_i; i < r_i; i++)
                {
                    vtkm::Id index = k * m_xdim * m_ydim + j * m_xdim + i;
                    // if (workIndex == 0)
                    //{
                    //     std::cout << index << std::endl;
                    // }
                    //  access the global array
                    OutputType originalvalue = inPointFieldPortal.Get(index);
                    boxMin = vtkm::Min(boxMin, originalvalue);
                    boxMax = vtkm::Max(boxMax, originalvalue);
                }
            }
        }

        minValue = boxMin;
        maxValue = boxMax;
    }

private:
    vtkm::IdComponent m_blockSize;
    vtkm::Id m_xdim;
    vtkm::Id m_ydim;
    vtkm::Id m_zdim;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN3D_h
