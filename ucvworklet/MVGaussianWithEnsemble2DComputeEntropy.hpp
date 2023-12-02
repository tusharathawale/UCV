#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_EL_COMPUTE_ENTROPY_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_EL_COMPUTE_ENTROPY_h

#include <vtkm/worklet/WorkletMapTopology.h>

class MVGaussianWithEnsemble2DComputeEntropy : public vtkm::worklet::WorkletMapField
{
public:
    MVGaussianWithEnsemble2DComputeEntropy(int numSamples) : m_numSamples(numSamples){};

    using ControlSignature = void(FieldInOut, FieldInOut, FieldInOut, WholeArrayIn);

    using ExecutionSignature = void(_1, _2, _3, _4, WorkIndex);

    // the first parameter is binded with the worklet
    using InputDomain = _1;

    template <typename CrossProbType, typename NonZeroProbType, typename EntropyType, typename CasesArrayType>
    VTKM_EXEC void operator()(
        CrossProbType &crossProb,
        NonZeroProbType &numNonzeroProb,
        EntropyType &entropy,
        CasesArrayType casesArray,
        vtkm::Id workIndex) const
    {
        // compute the start and end position for associated segment in the caseArray
        // the start index is the cellid*m_numSamples
        // the upper limit for the thread is the workIndex
        int startIndex = workIndex * m_numSamples;
        int endIndex = startIndex + m_numSamples - 1;
        // populate associated values into the probHist
        vtkm::Vec<vtkm::FloatDefault, 16> probHistogram(0);

        // printf("debug workindex %d startIndex %d endIndex %d\n",workIndex,startIndex,endIndex);

        // compute crossProb from probHist
        for (int index = startIndex; index <= endIndex; index++)
        {
            // it provides a readportal for whole array in
            vtkm::UInt8 caseValue = casesArray.Get(index);
            probHistogram[caseValue] = probHistogram[caseValue] + 1.0;
        }
        for (int i = 0; i < 16; i++)
        {
            probHistogram[i] = (probHistogram[i] / (1.0 * m_numSamples));
        }

        // cross probability
        // outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);
        crossProb = 1.0 - (probHistogram[0] + probHistogram[15]);

        vtkm::Id nonzeroProb = 0;
        vtkm::FloatDefault entropyValue = 0;
        vtkm::FloatDefault templog = 0;
        // compute number of nonzero cases
        // compute entropy
        for (int i = 0; i < 16; i++)
        {
            if (probHistogram[i] > 0.0001)
            {
                nonzeroProb++;
                templog = vtkm::Log2(probHistogram[i]);
                // if (i != 0 && i != totalNumCases - 1)
                //{
                //     totalnonzeroProb += probHistogram[i];
                // }
            }
            // do not update entropy if the pro is zero
            entropyValue = entropyValue + (-probHistogram[i]) * templog;
        }

        numNonzeroProb = nonzeroProb;
        entropy = entropyValue;
    }

private:
    int m_numSamples;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h
