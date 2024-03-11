#ifndef UCV_EXTRACTING_MIN_MAX_BYERR_h
#define UCV_EXTRACTING_MIN_MAX_BYERR_h

#include <vtkm/worklet/WorkletMapField.h>
#include <math.h>
#include <float.h>

struct ExtractMinMaxByErr : public vtkm::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldOut, FieldOut);
    using ExecutionSignature = void(_1, _2, _3, WorkIndex);
    using InputDomain = _1;
    ExtractMinMaxByErr(vtkm::Float64 err) : Error(err){};
    template <typename OriginalValuesType, typename OutputType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputType &minValue, OutputType &maxValue, vtkm::Id WorkIndex) const
    {

        minValue = vtkm::Float64(originalValues * 1.0) - this->Error/2.0;
        maxValue = vtkm::Float64(originalValues * 1.0) + this->Error/2.0;
        //if (WorkIndex >=0 && WorkIndex <=10)
        //{
        //    printf("err is %lf %lf %lf\n",this->Error, minValue, vtkm::Float64(originalValues * 1.0) - 2.0 * this->Error);
        //    printf("debug %.6f %.6f %.6f\n", originalValues, minValue, maxValue);
        //}
    }

    vtkm::Float64 Error = 0.0;
};

#endif // UCV_EXTRACTING_MIN_MAX_h