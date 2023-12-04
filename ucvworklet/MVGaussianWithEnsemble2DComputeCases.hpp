#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_EL_COMPUTE_CASES_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_EL_COMPUTE_CASES_h

#include <vtkm/worklet/WorkletMapTopology.h>

class MVGaussianWithEnsemble2DComputeCases : public vtkm::worklet::WorkletMapField
{
public:
    MVGaussianWithEnsemble2DComputeCases(int numCells, int numSamples, double isoValue) : m_numCells(numCells), m_numSamples(numSamples), m_isovalue(isoValue){};

    using ControlSignature = void(FieldInOut, WholeArrayIn, WholeArrayIn, WholeArrayIn);

    using ExecutionSignature = void(_1, _2, _3, _4, WorkIndex);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename CaseValueType, typename MatrixArrayType, typename MeanArrayType, typename SampleArrayType>
    VTKM_EXEC void operator()(
        CaseValueType &caseValue,
        MatrixArrayType matrixArray,
        MeanArrayType meanArray,
        SampleArrayType sampleArray,
        vtkm::Id workIndex) const
    {
        // if (workIndex == 0)
        // {
        //     printf("matrix array len %d\n", matrixArray.GetNumberOfValues());
        //     printf("mean array len %d\n", meanArray.GetNumberOfValues());
        //     printf("sample array len %d\n", sampleArray.GetNumberOfValues());
        // }

        // select matrix mean sample according to the workIndex
        int cellIndex = workIndex / m_numSamples;
        int sampleIndex = workIndex % m_numSamples;

        //printf("debug workIndex %d cellIndex %d sampleIndex %d\n",workIndex,cellIndex,sampleIndex);

        // compute A*S+M and decide the cases
        vtkm::Matrix<vtkm::FloatDefault, 4, 4> A = matrixArray.Get(cellIndex);
        vtkm::Vec<vtkm::FloatDefault, 4> M = meanArray.Get(cellIndex);
        vtkm::Vec<vtkm::FloatDefault, 4> S = sampleArray.Get(sampleIndex);
        vtkm::Vec<vtkm::FloatDefault, 4> result(0);

        // TODO, filter out cellMean and cellMax
        // need another array to show if the current cell need to do sampling

        for (uint i = 0; i < 4; i++)
        {
            for (uint j = 0; j < 4; j++)
            {
                result[i] += A[i][j] * S[j];
            }
            result[i] = result[i] + M[i];
        }
        
        //printf("worklet index %d ASM result %f %f %f %f\n",workIndex,result[0],result[1],result[2],result[3]);

        // TODO, simplify previous operation, A*M can be merged together in previous worklet
        // then we only need one vector here, just + M
        // pay attention to set init value as 0
        uint tempCaseValue = 0;
        for (uint i = 0; i < 4; i++)
        {
            // setting associated position to 1 if iso larger then specific cases
            if (m_isovalue >= result[i])
            {
                tempCaseValue = (1 << i) | tempCaseValue;
            }
        }
        caseValue=tempCaseValue;
        //printf("worklet index %d case value %d\n",workIndex,caseValue);
    }

private:
    int m_numCells;
    int m_numSamples;
    double m_isovalue;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h
