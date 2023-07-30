#ifndef UCV_ENTROPY_KDE_h
#define UCV_ENTROPY_KDE_h

#include <vtkm/worklet/WorkletMapField.h>
#include <cmath>


class EntropyKDE : public vtkm::worklet::WorkletMapField
{
public:
    EntropyKDE(double isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(FieldIn,
                                  FieldOut,
                                  FieldOut);

    using ExecutionSignature = void(_1, _2, _3);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble,
              typename OutPointFieldType>
    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutPointFieldType &outProbPostive,
        OutPointFieldType &outProbNegative) const
    {
        //for each vertex how to process the case where there are multiple variables
        vtkm::IdComponent numPoints = inPointFieldVecEnsemble.GetNumberOfComponents();
        std::cout<<"the number of points associated with each vertex is "<< numPoints << std::endl;
    }

private:
    double m_isovalue;
};

#endif // UCV_ENTROPY_KDE_h
