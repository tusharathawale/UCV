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

#endif // UCV_COMPUTE_DIFF_SUM_h