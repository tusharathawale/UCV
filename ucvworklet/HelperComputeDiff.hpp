#ifndef UCV_HELPER_COMPUTE_DIFF_h
#define UCV_HELPER_COMPUTE_DIFF_h

#include <vtkm/worklet/WorkletMapField.h>
#include <cmath>

class HelperComputeDiff : public vtkm::worklet::WorkletMapField
{
public:

    HelperComputeDiff(){};
    using ControlSignature = void(FieldIn,
                                  FieldIn,
                                  FieldOut);

    using ExecutionSignature = void(_1, _2, _3);

    // InPointFieldType should be a vector
    template <typename InType,
              typename OutType>
    VTKM_EXEC void operator()(
        const InType &field1,
        const InType &field2,
        OutType &diff) const
    {
        diff = vtkm::Abs(field1 - field2);
    }
};

#endif // UCV_HELPER_COMPUTE_DIFF_h
