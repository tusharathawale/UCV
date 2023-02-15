#ifndef UCV_EXTRACTING_MEAD_STD_h
#define UCV_EXTRACTING_MEAD_STD_h

#include <vtkm/worklet/WorkletReduceByKey.h>
#include <cmath>

struct ExtractingMeanStdev : public vtkm::worklet::WorkletReduceByKey
{
    using ControlSignature = void(KeysIn, ValuesIn, ReducedValuesOut, ReducedValuesOut);
    using ExecutionSignature = void(_2, _3, _4);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputType &meanValue, OutputType &stdevValue) const
    {
        vtkm::IdComponent NumComponents = originalValues.GetNumberOfComponents();

        // refer to https://www.strchr.com/standard_deviation_in_one_pass
        meanValue = originalValues[0];
        stdevValue = originalValues[0] * originalValues[0];
        for (vtkm::IdComponent index = 1; index < NumComponents; index++)
        {
            meanValue += originalValues[index];
            stdevValue += originalValues[index] * originalValues[index];
        }

        meanValue *= 1.0 / NumComponents;
        stdevValue *= 1.0 / NumComponents;
        stdevValue -= meanValue * meanValue;
        stdevValue = vtkm::Sqrt(stdevValue);
    }
};

#endif // UCV_EXTRACTING_MEAD_STD_h
