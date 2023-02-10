#ifndef UCV_ENTROPY_INDEPEDENT_GAUSSIAN_h
#define UCV_ENTROPY_INDEPEDENT_GAUSSIAN_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>
// compute the entropy and other assocaited uncertainty values *per cell*
class EntropyIndependentGaussian : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    EntropyIndependentGaussian(double isovalue)
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
    template <typename InPointFieldMeanType, typename InPointFieldStdevType, typename OutCellFieldType1, typename OutCellFieldType2, typename OutCellFieldType3>

    VTKM_EXEC void operator()(
        const InPointFieldMeanType &inPointFieldVecMean,
        const InPointFieldStdevType &inPointFieldVecStdev,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecMean.GetNumberOfComponents();
        // there are 8 points for each cell
        if (numPoints != 8)
        {
            printf("this is the 3d version for 8 vertecies\n");
            return;
        }
        vtkm::FloatDefault allPositiveProb = 1.0;
        vtkm::FloatDefault allNegativeProb = 1.0;
        vtkm::FloatDefault allCrossProb = 0.0;

        vtkm::FloatDefault positiveProb = 0.0;
        vtkm::FloatDefault negativeProb = 0.0;

        vtkm::Vec<vtkm::Vec2f, 8> ProbList;

        // there are 2^n total cases
        // int totalNumCases = static_cast<int>(vtkm::Pow(2.0, static_cast<vtkm::FloatDefault>(numPoints)));
        int totalNumCases = 256;
        vtkm::Vec<vtkm::FloatDefault, 256> probHistogram;

        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            vtkm::FloatDefault mean = inPointFieldVecMean[pointIndex];
            vtkm::FloatDefault stdev = inPointFieldVecStdev[pointIndex];

            // assuming we use the indepedent gaussian distribution
            // this is the error function to compute Pr[X<=L(m_iso)] for gaussian distribution
            negativeProb = 0.5 * (1 + std::erf((m_isovalue - mean) / (std::sqrt(2) * stdev)));
            positiveProb = 1.0 - negativeProb;

            allPositiveProb *= positiveProb;
            allNegativeProb *= negativeProb;

            ProbList[pointIndex][0] = negativeProb;
            ProbList[pointIndex][1] = positiveProb;
        }

        allCrossProb = 1 - allPositiveProb - allNegativeProb;
        outCellFieldCProb = allCrossProb;

        // TODO, use the number of vertesies as another parameter
        // traverseRec(1.0, 0, 0, numPoints, ProbList, probHistogram);
        traverseBit(ProbList, probHistogram);
        // go through each option in the case table
        // 1 is positive 0 is negative
        // extracting the entropy or other values based on probHistogram

        vtkm::FloatDefault entropyValue = 0;
        vtkm::Id nonzeroCases = 0;
        vtkm::FloatDefault templog = 0;

        // use this to check the reuslts as needed
        // vtkm::FloatDefault totalnonzeroProb = 0;

        for (int i = 0; i < totalNumCases; i++)
        {
            templog = 0;
            if (probHistogram[i] > 0.00001)
            {
                nonzeroCases++;
                templog = vtkm::Log2(probHistogram[i]);
                //if (i != 0 && i != totalNumCases - 1)
                //{
                //    totalnonzeroProb += probHistogram[i];
                //}
            }
            entropyValue = entropyValue + (-probHistogram[i]) * templog;
        }
        outCellFieldNumNonzeroProb = nonzeroCases;
        outCellFieldEntropy = entropyValue;

        //if (allCrossProb != 0 || totalnonzeroProb != 0)
        //{
            // this is for correctness checking
            //if (fabs(allCrossProb - totalnonzeroProb) > 0.001)
            //{
                //std::cout << "bad value " << allCrossProb << " " << totalnonzeroProb << std::endl;
            //}
        //}

    }

    VTKM_EXEC inline void traverseBit(vtkm::Vec<vtkm::Vec2f, 8> &ProbList,
                                      vtkm::Vec<vtkm::FloatDefault, 256> &probHistogram) const
    {

        // go through each option in the case table
        // 1 is positive 0 is negative
        for (uint i = 0; i < 256; i++)
        {
            vtkm::FloatDefault currProb = 1.0;
            for (uint j = 0; j < 8; j++)
            {
                if (i & (1 << j))
                {
                    // positive case
                    currProb = currProb * ProbList[j][1];
                }
                else
                {
                    // negative case
                    currProb = currProb * ProbList[j][0];
                }
            }
            probHistogram[i] = currProb;
        }
    }

    // using recursive call to go through all possibilities
    void traverseRec(vtkm::FloatDefault currentProb, int depth, int id, const int numPoints,
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

        traverseRec(nextPosProb, depth + 1, 1 + (id << 1), numPoints, ProbList, probHistogram);
        traverseRec(nextNegProb, depth + 1, id << 1, numPoints, ProbList, probHistogram);
        return;
    }

private:
    double m_isovalue;
};

#endif // UCV_ENTROPY_INDEPEDENT_GAUSSIAN_h