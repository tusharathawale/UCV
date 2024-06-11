
#ifndef UCV_COMPUTE_DIFF_h
#define UCV_COMPUTE_DIFF_h

#include <vtkm/worklet/WorkletMapField.h>
#include <cmath>

// go through each vertexies and compute associated mean and stdec
struct ComputeDiff: public vtkm::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldIn, FieldOut);
    using ExecutionSignature = void(_1, _2, _3, WorkIndex);
    template <typename TypeIn, typename TypeOut>
    VTKM_EXEC void operator()(
        const TypeIn &input1, const TypeIn &input2, TypeOut &diff, vtkm::Id workIndex) const
    {
        //not sure if this is essential
        diff = vtkm::Abs(input1-input2);
        if(diff > 0.5){
            ("index %ld\n",workIndex);
        }
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