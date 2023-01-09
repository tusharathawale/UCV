#ifndef UCV_EXTRACTING_MIN_MAX_h
#define UCV_EXTRACTING_MIN_MAX_h

#include <vtkm/worklet/WorkletReduceByKey.h>

struct ExtractingMinMax : public vtkm::worklet::WorkletReduceByKey
{
    using ControlSignature = void(KeysIn, ValuesIn, ReducedValuesOut, ReducedValuesOut);
    using ExecutionSignature = void(_2, _3, _4);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputType &minValue, OutputType &maxValue) const
    {
        vtkm::FloatDefault boxMin = DBL_MAX;
        vtkm::FloatDefault boxMax = 0;

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