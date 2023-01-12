#ifndef UCV_EXTRACTING_MEAN_RAW_h
#define UCV_EXTRACTING_MEAN_RAW_h

#include <vtkm/worklet/WorkletReduceByKey.h>

struct ExtractingMeanRaw : public vtkm::worklet::WorkletReduceByKey
{
    using ControlSignature = void(KeysIn, ValuesIn, ReducedValuesOut, ReducedValuesOut);
    using ExecutionSignature = void(_2, _3, _4);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputScalarType, typename OutputVecType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputScalarType &meanValue, OutputVecType &rawVec) const
    {
        vtkm::FloatDefault boxSum = 0;
        vtkm::Id NumComponents = originalValues.GetNumberOfComponents();

        for (vtkm::IdComponent index = 0;
             index < originalValues.GetNumberOfComponents(); index++)
        {
            boxSum = boxSum + static_cast<vtkm::FloatDefault>(originalValues[index]);
            rawVec[index] = originalValues[index];
        }

        meanValue = boxSum / (1.0 * NumComponents);
    }
};

#endif // UCV_EXTRACTING_MEAN_RAW_h