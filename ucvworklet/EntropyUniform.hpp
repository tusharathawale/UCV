#ifndef UCV_ENTROPY_UNIFORM_h
#define UCV_ENTROPY_UNIFORM_h

#include <vtkm/worklet/WorkletMapTopology.h>
class EntropyUniform : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    EntropyUniform(double isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4, _5, _6);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldMinType, typename InPointFieldMaxType, typename OutCellFieldType1, typename OutCellFieldType2, typename OutCellFieldType3>

    VTKM_EXEC void operator()(
        const InPointFieldMinType &inPointFieldVecMin,
        const InPointFieldMaxType &inPointFieldVecMax,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecMin.GetNumberOfComponents();
        // there are 8 points for each cell

        if (numPoints != 8)
        {
            printf("this is the 3d version for 8 vertecies\n");
            return;
        }

        vtkm::FloatDefault allPositiveProb = 1.0;
        vtkm::FloatDefault allNegativeProb = 1.0;
        vtkm::FloatDefault allCrossProb = 0.0;

        vtkm::FloatDefault positiveProb;
        vtkm::FloatDefault negativeProb;

        // position 0 is negative
        // position 1 is positive
        vtkm::Vec<vtkm::Vec2f, 8> ProbList;

        // there are 2^n total cases
        // int totalNumCases = static_cast<int>(vtkm::Pow(2.0, static_cast<vtkm::FloatDefault>(numPoints)));
        int totalNumCases = 256;
        // std::vector<vtkm::FloatDefault> probHistogram;
        // probHistogram.resize(totalNumCases);
        //  std::cout << "debug totalNumCases " << totalNumCases << std::endl;

        vtkm::Vec<vtkm::FloatDefault, 256> probHistogram;

        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            vtkm::FloatDefault minV = inPointFieldVecMin[pointIndex];
            vtkm::FloatDefault maxV = inPointFieldVecMax[pointIndex];

            if (this->m_isovalue <= minV)
            {
                positiveProb = 1.0;
                negativeProb = 0.0;
            }
            else if (this->m_isovalue >= maxV)
            {
                positiveProb = 0.0;
                negativeProb = 1.0;
            }
            else
            {
                // assuming we use the uniform distribution
                positiveProb = (maxV - this->m_isovalue) / (maxV - minV);
                negativeProb = 1.0 - positiveProb;
            }

            // positiveProbList[pointIndex] = positiveProb;
            // negativeProbList[pointIndex] = negativeProb;
            allNegativeProb *= negativeProb;
            allPositiveProb *= positiveProb;

            ProbList[pointIndex][0] = negativeProb;
            ProbList[pointIndex][1] = positiveProb;
        }

        allCrossProb = 1 - allPositiveProb - allNegativeProb;
        outCellFieldCProb = allCrossProb;

        // printf("debug cuda, ok allCrossProb\n");

        // TODO, use the number of vertesies as another parameter
        // there is recursion call and the nvlink might give warning
        // such as the stack size can not be determined statically
        // traverse(1.0, 0, 0, numPoints, ProbList, probHistogram);

        traverseBit(ProbList, probHistogram);

        // extracting the entropy or other values based on probHistogram
        vtkm::FloatDefault entropyValue = 0;
        vtkm::Id nonzeroCases = 0;
        vtkm::FloatDefault templog = 0;

        // printf("debug cuda, ok probHistogram\n");

        // test
        vtkm::FloatDefault totalnonzeroProb = 0;

        for (int i = 0; i < totalNumCases; i++)
        {
            templog = 0;
            if (probHistogram[i] > 0.00001)
            {
                nonzeroCases++;
                templog = vtkm::Log2(probHistogram[i]);
                if (i != 0 && i != totalNumCases - 1)
                {
                    totalnonzeroProb += probHistogram[i];
                }
            }
            entropyValue = entropyValue + (-probHistogram[i]) * templog;
        }

        outCellFieldNumNonzeroProb = nonzeroCases;
        outCellFieldEntropy = entropyValue;

        if (allCrossProb != 0 || totalnonzeroProb != 0)
        {
            if (fabs(allCrossProb - totalnonzeroProb) > 0.001)
            {
                //std::cout << "bad value " << allCrossProb << " " << totalnonzeroProb << std::endl;
            }
        }
        // printf("debug cuda, ok entropy\n");
    }

    VTKM_EXEC inline void traverseBit(vtkm::Vec<vtkm::Vec2f, 8> &ProbList,
                                      vtkm::Vec<vtkm::FloatDefault, 256> &probHistogram) const
    {

        // go through each option in the case table
        // 1 is positive 0 is negative
        for (vtkm::UInt16 i = 0; i < 256; i++)
        {
            vtkm::FloatDefault currProb = 1.0;
            for (vtkm::UInt8 j = 0; j < 8; j++)
            {
                if (i & (1 << j))
                {
                    // positive case
                    currProb *= ProbList[j][1];
                }
                else
                {
                    // negative case
                    currProb *= ProbList[j][0];
                }
            }
            probHistogram[i] = currProb;
        }
    }

    // using recursive call to go through all possibilities
    // there are some cuda memory issue with recursive call here
    VTKM_EXEC inline void traverse(vtkm::FloatDefault currentProb, int depth, int id, const int numPoints,
                                   vtkm::Vec<vtkm::Vec2f, 8> &ProbList,
                                   vtkm::Vec<vtkm::FloatDefault, 256> &probHistogram) const
    {
        // TODO, make this as a private variable
        // how to set it as a private variable of the worklet
        if (depth == numPoints)
        {
            // if (id > 256) how to set this as a worklet parameter?
            //{
            //     throw std::runtime_error("id is supposed to be 0 to 255");
            // }
            probHistogram[id] = currentProb;
            return;
        }
        // two branches for current node
        vtkm::FloatDefault nextPosProb = currentProb * ProbList[depth][1];
        vtkm::FloatDefault nextNegProb = currentProb * ProbList[depth][0];

        traverse(nextPosProb, depth + 1, 1 + (id << 1), numPoints, ProbList, probHistogram);
        traverse(nextNegProb, depth + 1, id << 1, numPoints, ProbList, probHistogram);
        return;
    }

private:
    double m_isovalue;
};

#endif // UCV_ENTROPY_UNIFORM_h