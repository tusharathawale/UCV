#ifndef UCV_COMPUTE_DIFF_SUM_h
#define UCV_COMPUTE_DIFF_SUM_h

#include <vtkm/worklet/WorkletMapField.h>
#include <cmath>

// go through each vertexies and compute associated mean and stdec
struct ComputeDiffSum : public vtkm::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldIn, FieldOut);
    using ExecutionSignature = void(_1, _2, _3);
    template <typename TypeIn, typename TypeOut>
    VTKM_EXEC void operator()(
        const TypeIn &currEntropy, const TypeIn &allEntropy, TypeOut &diff) const
    {
        diff = vtkm::Abs(currEntropy-allEntropy);
        return;
    }
};

struct ComputeDiffSquare : public vtkm::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldIn, FieldOut);
    using ExecutionSignature = void(_1, _2, _3);
    template <typename TypeIn, typename TypeOut>
    VTKM_EXEC void operator()(
        const TypeIn &v1, const TypeIn &v2, TypeOut &diff) const
    {
        diff = (v1-v2)*(v1-v2);
        return;
    }
};

#endif // UCV_COMPUTE_DIFF_SUM_h