#ifndef UCV_EXTRACTING_MIN_MAX_FROM_MEAN_DEV_h
#define UCV_EXTRACTING_MIN_MAX_FROM_MEAN_DEV_h

#include <vtkm/worklet/WorkletMapField.h>

struct ExtractingMinMaxFromMeanDev : public vtkm::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldIn, FieldOut, FieldOut);
    using ExecutionSignature = void(_1, _2, _3, _4);
    using InputDomain = _1;
    template <typename InputType, typename OutputType>
    VTKM_EXEC void operator()(
        const InputType &mean, const InputType &dev, OutputType &minValue, OutputType &maxValue) const
    {
        minValue = mean - 1*dev;
        maxValue = mean + 1*dev;
    }
};

#endif // UCV_EXTRACTING_MIN_MAX_FROM_MEAN_DEV_h