#ifndef UCV_CRITICAL_POINT_h
#define UCV_CRITICAL_POINT_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>

struct CriticalPointWorklet : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    CriticalPointWorklet(){};

    using ControlSignature = void(CellSetIn, FieldIn);

    using ExecutionSignature = void(_2, Boundary, WorkIndex);

    template <typename InPointField>
    VTKM_EXEC void operator()(const InPointField &input,
                              const vtkm::exec::BoundaryState &boundary,
                              vtkm::Id WorkIndex) const
    {

        auto minIndices = boundary.MinNeighborIndices(this->m_neighborhoodSize);
        auto maxIndices = boundary.MaxNeighborIndices(this->m_neighborhoodSize);

        if (WorkIndex == 0)
        {
            // debug
            printf("workIndex is %d\n", WorkIndex);
            printf("min index %d %d %d\n", minIndices[0], minIndices[1], minIndices[2]);
            printf("max index %d %d %d\n", maxIndices[0], maxIndices[1], maxIndices[2]);
        }
    }

private:
    int m_neighborhoodSize = 1;
};

#endif // UCV_CRITICAL_POINT_h