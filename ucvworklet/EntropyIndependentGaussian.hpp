#ifndef UCV_ENTROPY_INDEPEDENT_GAUSSIAN_h
#define UCV_ENTROPY_INDEPEDENT_GAUSSIAN_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>
// TODO, update this into indepedent gaussian version
class EntropyIndependentGaussian : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    EntropyIndependentGaussian(int isovalue)
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

        vtkm::FloatDefault allPositiveProb = 1.0;
        vtkm::FloatDefault allNegativeProb = 1.0;
        vtkm::FloatDefault allCrossProb = 0.0;

        vtkm::FloatDefault positiveProb;
        vtkm::FloatDefault negativeProb;

        std::vector<vtkm::FloatDefault> positiveProbList;
        std::vector<vtkm::FloatDefault> negativeProbList;

        positiveProbList.resize(numPoints);
        negativeProbList.resize(numPoints);

        // there are 2^n total cases
        int totalNumCases = static_cast<int>(vtkm::Pow(2.0, static_cast<vtkm::FloatDefault>(numPoints)));
        std::vector<vtkm::FloatDefault> probHistogram;

        probHistogram.resize(totalNumCases);


        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            // TODO, computing mean and std based on the data in the neigoborhood
            vtkm::FloatDefault mean = inPointFieldVecMean[pointIndex];
            vtkm::FloatDefault stdev = inPointFieldVecStdev[pointIndex];
            
            // is this necessary for using min and max here?
            //if (this->m_isovalue <= minV)
            //{
            //    positiveProb = 1.0;
            //    negativeProb = 0.0;
            //}
            //else if (this->m_isovalue >= maxV)
            //{
            //    positiveProb = 0.0;
            //    negativeProb = 1.0;
            //}
            //else
            //{
                // assuming we use the indepedent gaussian distribution
                negativeProb = 0.5*(1 + std::erf((m_isovalue - mean)/(std::sqrt(2)*stdev)));
                positiveProb = 1.0 - negativeProb;
            //}

            positiveProbList[pointIndex] = positiveProb;
            negativeProbList[pointIndex] = negativeProb;

            allPositiveProb *= positiveProb;
            allNegativeProb *= negativeProb;
        }

        allCrossProb = 1 - allPositiveProb - allNegativeProb;
        outCellFieldCProb = allCrossProb;

        // TODO, use the number of vertesies as another parameter
        traverse(1.0, 0, 0, numPoints, positiveProbList, negativeProbList, probHistogram);

        // extracting the entropy or other values based on probHistogram

        vtkm::FloatDefault entropyValue = 0;
        vtkm::Id nonzeroCases = 0;
        vtkm::FloatDefault templog = 0;

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
            //std::cout << "test " << allCrossProb << " " << totalnonzeroProb << std::endl;
        }
    }

    // using recursive call to go through all possibilities
    void traverse(vtkm::FloatDefault currentProb, int depth, int id, const int numPoints,
                  std::vector<vtkm::FloatDefault> &positiveProbList,
                  std::vector<vtkm::FloatDefault> &negativeProbList,
                  std::vector<vtkm::FloatDefault> &probHistogram) const
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
        vtkm::FloatDefault nextPosProb = currentProb * positiveProbList[depth];
        vtkm::FloatDefault nextNegProb = currentProb * negativeProbList[depth];

        traverse(nextPosProb, depth + 1, 1 + (id << 1), numPoints, positiveProbList, negativeProbList, probHistogram);
        traverse(nextNegProb, depth + 1, id << 1, numPoints, positiveProbList, negativeProbList, probHistogram);
        return;
    }

private:
    int m_isovalue;
};


#endif //UCV_ENTROPY_INDEPEDENT_GAUSSIAN_h