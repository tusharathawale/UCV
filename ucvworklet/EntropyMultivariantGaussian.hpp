#ifndef UCV_ENTROPY_MULTIVARIANT_GAUSSIAN_h
#define UCV_ENTROPY_MULTIVARIANT_GAUSSIAN_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include "eigenmvn.h"

class EntropyMultivariantGaussian : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    EntropyUniform(double isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldMeanType, typename InPointFieldStdevType, typename OutCellFieldType1>

    VTKM_EXEC void operator()(
        const InPointFieldMeanType &inPointFieldVecMean,
        const InPointFieldStdevType &inPointFieldVecStdev,
        OutCellFieldType1 &outCellFieldCProb) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecMin.GetNumberOfComponents();

        if (numPoints != 4)
        {
            throw std::runtime_error("we only support 2d case for the multivariant gaussian distribution");
        }

        vtkm::FloatDefault allPositiveProb = 1.0;
        vtkm::FloatDefault allNegativeProb = 1.0;
        vtkm::FloatDefault allCrossProb = 0.0;

        // generate mean and cov matrix
        Eigen::Vector4d mean4by4;

        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            vtkm::FloatDefault mean = inPointFieldVecMean[pointIndex];
            vtkm::FloatDefault stdev = inPointFieldVecStdev[pointIndex];
            //TODO
        }

        allCrossProb = 1 - allPositiveProb - allNegativeProb;
        outCellFieldCProb = allCrossProb;
    }

private:
    double m_isovalue;
};

#endif // UCV_ENTROPY_MULTIVARIANT_GAUSSIAN_h