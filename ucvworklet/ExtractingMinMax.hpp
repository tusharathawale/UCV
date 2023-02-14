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
      minValue = originalValues[0];
      maxValue = originalValues[0];

      for (vtkm::IdComponent index = 1;
           index < originalValues.GetNumberOfComponents(); index++)
      {
        auto originalValue = originalValues[index];
        // This iteration is to support ArrayHandleRecombineVec from CastAndCallWithExtractedArray
        for (vtkm::IdComponent cIndex = 0; cIndex < originalValue.GetNumberOfComponents(); ++cIndex)
        {
          minValue[cIndex] =
              (originalValues[index][cIndex] < minValue[cIndex]) ? originalValues[index][cIndex] : minValue[cIndex];
          maxValue[cIndex] =
              (originalValues[index][cIndex] > maxValue[cIndex]) ? originalValues[index][cIndex] : maxValue[cIndex];
        }
      }
    }
};


#endif //UCV_EXTRACTING_MIN_MAX_h
