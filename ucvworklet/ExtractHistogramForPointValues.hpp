#ifndef UCV_EXTRACT_HISTOGRAM_OFPOINT_VALUES_h
#define UCV_EXTRACT_HISTOGRAM_OFPOINT_VALUES_h

#include <vtkm/worklet/WorkletMapField.h>
#include <math.h>
#include <float.h>

#ifdef USE_LOG
#define LOG(x) x
#else
#define LOG(x)
#endif

struct ExtractHistogramForPointValues : public vtkm::worklet::WorkletMapField
{
    ExtractHistogramForPointValues(vtkm::Id numBins) : m_numBins(numBins){};

    using ControlSignature = void(FieldIn, FieldOut, FieldOut);
    using ExecutionSignature = void(_1, _2, _3, WorkIndex);
    using InputDomain = _1;
    template <typename OriginalValuesType, typename OutputType>
    VTKM_EXEC void operator()(
        const OriginalValuesType &originalValues, OutputType &histDensity, OutputType &histEdges, vtkm::Id WorkIndex) const
    {

        // check the histogram size
        if (histDensity.GetNumberOfComponents() != this->m_numBins)
        {
            // the number of values in histogram vec should be same with num bins
            printf("Error, the component in histogram vec is supposed to same with numBins");
            return;
        }

        // extract min and max of input and output
        vtkm::FloatDefault min = DBL_MAX;
        vtkm::FloatDefault max = vtkm::NegativeInfinity64();

        for (vtkm::IdComponent index = 0;
             index < originalValues.GetNumberOfComponents(); index++)
        {
            vtkm::FloatDefault originalvalue = originalValues[index];
            min = vtkm::Min(min, originalvalue);
            max = vtkm::Max(max, originalvalue);
        }

        LOG(printf("debug index %d, first two values %lf, %lf\n", WorkIndex, static_cast<double>(originalValues[0]), static_cast<double>(originalValues[1]));)

        // compute the bin length
        vtkm::FloatDefault binSize = (max - min) / (1.0 * this->m_numBins);

        // compute the hist edges
        for (vtkm::IdComponent index = 0;
             index < histEdges.GetNumberOfComponents(); index++)
        {
            histEdges[index] = min + index * binSize;
        }

        // init the histDensity as zero
        for (vtkm::IdComponent index = 0;
             index < this->m_numBins; index++)
        {
            histDensity[index] = 0;
        }

        // go through each ensemble values to compute the histDensity
        for (vtkm::IdComponent index = 0;
             index < originalValues.GetNumberOfComponents(); index++)
        {
            vtkm::FloatDefault originalvalue = originalValues[index];
            vtkm::Id binIndex = static_cast<vtkm::Id>((originalvalue - min) / binSize);
            // consider the largest element, it should be included into the last bin slot
            if (vtkm::Abs(originalvalue - max) < 0.0000001 && binIndex == this->m_numBins)
            {
                binIndex = binIndex - 1;
            }
            histDensity[binIndex] += 1;
        }

        // normalize the histDensity
        vtkm::FloatDefault histSum = 0;
        for (vtkm::IdComponent index = 0;
             index < this->m_numBins; index++)
        {
            histSum += histDensity[index];
        }

        for (vtkm::IdComponent index = 0;
             index < this->m_numBins; index++)
        {
            histDensity[index] = histDensity[index] / histSum;
        }
    }

    vtkm::Id m_numBins = 5;
};

#endif // UCV_EXTRACT_HISTOGRAM_OFPOINT_VALUES_h