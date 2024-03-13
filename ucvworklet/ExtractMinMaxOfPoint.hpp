#ifndef UCV_EXTRACTING_MIN_MAX_OFPOINT_h
#define UCV_EXTRACTING_MIN_MAX_OFPOINT_h

#include <vtkm/worklet/WorkletMapField.h>
#include <math.h>
#include <float.h>

struct ExtractMinMaxOfPoint : public vtkm::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldOut, FieldOut);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputType &minValue, OutputType &maxValue) const
    {
        vtkm::FloatDefault boxMin = DBL_MAX;
        //the data can be a negative number
        vtkm::FloatDefault boxMax = vtkm::NegativeInfinity64();

        for (vtkm::IdComponent index = 0;
             index < originalValues.GetNumberOfComponents(); index++)
        {
            vtkm::FloatDefault originalvalue = originalValues[index];
            boxMin = vtkm::Min(boxMin, originalvalue);
            boxMax = vtkm::Max(boxMax, originalvalue);
        }

        minValue = boxMin;
        maxValue = boxMax;
    }
};


#endif //UCV_EXTRACTING_MIN_MAX_h