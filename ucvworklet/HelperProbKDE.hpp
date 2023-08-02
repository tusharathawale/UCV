#ifndef UCV_Prob_KDE_h
#define UCV_Prob_KDE_h

#include <vtkm/worklet/WorkletMapField.h>
#include <cmath>
#include "./linalg/EasyLinalg/kde.h"

class HelperProbKDE : public vtkm::worklet::WorkletMapField
{
public:
    HelperProbKDE(double isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(FieldIn,
                                  FieldOut);

    using ExecutionSignature = void(_1, _2);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble,
              typename OutPointFieldType>
    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutPointFieldType &outProbPostive) const
    {
        // for each vertex how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecEnsemble.GetNumberOfComponents();
        // std::cout<<"the number of points associated with each vertex is "<< numPoints << std::endl;
        // the numPoints has 64 values
        
        // traverse the inputPointFieldIntoTheEasyLinalg and change it into the dedicated format
        EASYLINALG::Vec<float, 4 * 4 * 4> ensembles;
        for (vtkm::IdComponent pointIndex = 0; pointIndex < numPoints; ++pointIndex)
        {
            ensembles[pointIndex] = inPointFieldVecEnsemble[pointIndex];
        }

        // using this value to compute kde
        outProbPostive = EASYLINALG::KDECDF1D(ensembles, this->m_isovalue);
    }

private:
    float m_isovalue;
};

#endif // UCV_Prob_KDE_h
