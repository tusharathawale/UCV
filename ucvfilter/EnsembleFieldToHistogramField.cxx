#include <EnsembleFieldToHistogramField.h>

#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/Math.h>
#include <vtkm/Range.h>

#ifdef USE_LOG
#define LOG(x) x
#else
#define LOG(x)
#endif

namespace
{

struct ExtractHistogramForPointValues : public vtkm::worklet::WorkletMapField
{
  ExtractHistogramForPointValues(vtkm::Id numBins)
    : m_numBins(numBins){};

  using ControlSignature = void(FieldIn, FieldOut, FieldOut);
  using ExecutionSignature = void(_1, _2, _3, WorkIndex);
  using InputDomain = _1;
  template <typename OriginalValuesType, typename OutputType>
  VTKM_EXEC void operator()(const OriginalValuesType& originalValues,
                            OutputType& histDensity,
                            OutputType& histEdges,
                            vtkm::Id LOG(WorkIndex)) const
  {
    using ComponentType = typename OutputType::ComponentType;

    // check the histogram size
    if (histDensity.GetNumberOfComponents() != this->m_numBins)
    {
      // the number of values in histogram vec should be same with num bins
      printf("Error, the component in histogram vec is supposed to same with numBins");
      return;
    }

    // extract min and max of input and output
    vtkm::Range range;

    for (vtkm::IdComponent index = 0; index < originalValues.GetNumberOfComponents(); index++)
    {
      range.Include(originalValues[index]);
    }

    LOG(printf("debug index %d, first two values %lf, %lf\n",
               WorkIndex,
               static_cast<double>(originalValues[0]),
               static_cast<double>(originalValues[1]));)

    // compute the bin length
    vtkm::Float64 binSize = range.Length() / this->m_numBins;

    // compute the hist edges
    for (vtkm::IdComponent index = 0; index < histEdges.GetNumberOfComponents(); index++)
    {
      histEdges[index] = static_cast<ComponentType>(range.Min + index * binSize);
    }

    // init the histDensity as zero
    for (vtkm::IdComponent index = 0; index < this->m_numBins; index++)
    {
      histDensity[index] = 0;
    }

    // go through each ensemble values to compute the histDensity
    for (vtkm::IdComponent index = 0; index < originalValues.GetNumberOfComponents(); index++)
    {
      vtkm::FloatDefault originalvalue = originalValues[index];
      vtkm::Id binIndex = static_cast<vtkm::Id>((originalvalue - range.Min) / binSize);
      // consider the largest element, it should be included into the last bin slot
      if (binIndex >= this->m_numBins)
      {
        binIndex = binIndex - 1;
      }
      histDensity[binIndex] += 1;
    }

    // normalize the histDensity
    ComponentType histSum = 0;
    for (vtkm::IdComponent index = 0; index < this->m_numBins; index++)
    {
      histSum += histDensity[index];
    }

    ComponentType weight = 1 / histSum;
    for (vtkm::IdComponent index = 0; index < this->m_numBins; index++)
    {
      histDensity[index] = histDensity[index] * weight;
    }
  }

  vtkm::Id m_numBins = 5;
};

} // anonymous namespace

namespace vtkm
{
namespace filter
{
namespace uncertainty
{

VTKM_CONT EnsembleFieldToHistogramField::EnsembleFieldToHistogramField()
{
  this->SetHistogramDensityName("HistogramDensity");
}

VTKM_CONT vtkm::cont::DataSet EnsembleFieldToHistogramField::DoExecute(
  const vtkm::cont::DataSet& inData)
{
  vtkm::cont::Field ensembleField = this->GetFieldFromDataSet(inData);
  vtkm::cont::UnknownArrayHandle ensembleUnknown = ensembleField.GetData();
  vtkm::IdComponent ensembleComponents = ensembleUnknown.GetNumberOfComponentsFlat();
  if (ensembleComponents < 2)
  {
    throw vtkm::cont::ErrorBadType("Cannot compute probabilities with only one ensemble member.");
  }

  vtkm::IdComponent numBins = this->GetNumberOfBins();
  if (numBins < 2)
  {
    // Use Sturge's rule for histogram size.
    numBins =
      1 + static_cast<vtkm::IdComponent>(vtkm::Round(3.322 * vtkm::Log10(ensembleComponents)));
  }

  vtkm::cont::UnknownArrayHandle histogramDensityData;
  vtkm::cont::UnknownArrayHandle histogramEdgesData;

  auto resolveType = [&](auto ensemble)
  {
    using ComponentType = typename decltype(ensemble)::ValueType::ComponentType;

    vtkm::cont::ArrayHandleRuntimeVec<ComponentType> histogramDensity(numBins);
    vtkm::cont::ArrayHandleRuntimeVec<ComponentType> histogramEdges(numBins + 1);

    this->Invoke(
      ExtractHistogramForPointValues{ numBins }, ensemble, histogramDensity, histogramEdges);

    histogramDensityData = histogramDensity;
    histogramEdgesData = histogramEdges;
  };

  ensembleUnknown.CastAndCallWithExtractedArray(resolveType);

  vtkm::cont::DataSet outData = this->CreateResultField(
    inData, this->GetHistogramDensityName(), ensembleField.GetAssociation(), histogramDensityData);
  outData.AddField(
    { this->GetHistogramEdgesName(), ensembleField.GetAssociation(), histogramEdgesData });
  return outData;
}

}
}
} // namespace vtkm::filter::uncertainty
