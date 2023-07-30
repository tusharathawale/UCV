#ifndef UCV_EXTRACTING_RAW_h
#define UCV_EXTRACTING_RAW_h

#include <vtkm/worklet/WorkletReduceByKey.h>

struct ExtractingRaw : public vtkm::worklet::WorkletReduceByKey
{
    using ControlSignature = void(KeysIn, ValuesIn, ReducedValuesOut);
    using ExecutionSignature = void(_2, _3);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputVecType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputVecType &rawVec) const
    {
        vtkm::FloatDefault boxSum = 0;
        vtkm::Id NumComponents = originalValues.GetNumberOfComponents();

        for (vtkm::IdComponent index = 0;
             index < originalValues.GetNumberOfComponents(); index++)
        {
            rawVec[index] = originalValues[index];
        }
    }
};

#endif // UCV_EXTRACTING_RAW_h