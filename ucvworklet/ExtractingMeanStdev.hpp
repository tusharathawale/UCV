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
        vtkm::FloatDefault boxSum = 0;

        vtkm::IdComponent NumComponents = originalValues.GetNumberOfComponents();

        // refer to https://www.strchr.com/standard_deviation_in_one_pass
        for (vtkm::IdComponent index = 0;
             index < NumComponents; index++)
        {
            boxSum = boxSum + static_cast<vtkm::FloatDefault>(originalValues[index]);
        }

        meanValue = boxSum / (1.0 * (NumComponents));

        vtkm::FloatDefault diffSum = 0;

        for (vtkm::IdComponent index = 0;
             index < NumComponents; index++)
        {
            vtkm::FloatDefault diff = static_cast<vtkm::FloatDefault>(originalValues[index]) - static_cast<vtkm::FloatDefault>(meanValue);
            diffSum += diff * diff;
        }

        stdevValue = std::sqrt(diffSum / (1.0*NumComponents));
    }
};

#endif // UCV_EXTRACTING_MEAD_STD_h