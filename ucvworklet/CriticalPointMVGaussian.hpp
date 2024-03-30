#ifndef UCV_CRITICAL_POINT_MV_GAUSSIAN_h
#define UCV_CRITICAL_POINT_MV_GAUSSIAN_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>
#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
#include <thrust/device_vector.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/uniform_real_distribution.h>
#else
#include <random>
#endif

#ifdef USE_LOG
#define LOG(x) x
#else
#define LOG(x)
#endif

struct CriticalPointMVGaussian : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    CriticalPointMVGaussian(vtkm::Id samples) : m_NumSamples(samples){};

    using ControlSignature = void(CellSetIn, FieldInNeighborhood, FieldOut);

    using ExecutionSignature = void(_2, _3, Boundary);

    template <typename InPointField, typename OutPointField>
    VTKM_EXEC void operator()(const InPointField &inPointFieldVecEnsemble,
                              OutPointField &minProb,
                              const vtkm::exec::BoundaryState &boundary) const
    {
        // resluts is the coordinates of three dims
        auto minIndices = boundary.MinNeighborIndices(1);
        auto maxIndices = boundary.MaxNeighborIndices(1);

        // minIndices is supposed to be -1
        // maxIndices is supposed to be 1
        // if (WorkIndex == 0)
        //{
        // debug
        LOG(printf("workIndex is %d\n", WorkIndex));
        // printf("min index %d %d %d\n", minIndices[0], minIndices[1], minIndices[2]);
        // printf("max index %d %d %d\n", maxIndices[0], maxIndices[1], maxIndices[2]);

        // filter out the element in the boundry
        // if the element is at the boundry places its min prob is 0
        if ((maxIndices[0] - minIndices[0] < 2) || (maxIndices[1] - minIndices[1] < 2))
        {
            // if x and y is at the boundry, do not consider it
            minProb = 0;
            return;
        }

        using VecType = decltype(inPointFieldVecEnsemble.Get(0, 0, 0));

        //vtkm::FloatDefault m1 = minValue.Get(0, 0, 0);
        VecType ensArray1 = inPointFieldVecEnsemble.Get(0, 0, 0);

        print("get number of component %d\n", ensArray1.GetNumberOfComponents());

        // vtkm::FloatDefault m2 = minValue.Get(0, 1, 0);
        // vtkm::FloatDefault ensArray2 = inPointFieldVecEnsemble.Get(0, 1, 0);

        // vtkm::FloatDefault m3 = minValue.Get(0, -1, 0);
        // vtkm::FloatDefault ensArray3 = inPointFieldVecEnsemble.Get(0, -1, 0);

        // vtkm::FloatDefault m4 = minValue.Get(1, 0, 0);
        // vtkm::FloatDefault ensArray4 = inPointFieldVecEnsemble.Get(1, 0, 0);

        // vtkm::FloatDefault m5 = minValue.Get(-1, 0, 0);
        // vtkm::FloatDefault ensArray5 = inPointFieldVecEnsemble.Get(-1, 0, 0);

        // if (abs(m2-0.0)<0.0000001 || abs(m3-0.0)<0.0000001 || abs(m4-0.0)<0.0000001 || abs(m5-0.0)<0.0000001 ){
        //     minProb = 0;
        //     return;
        // }

        int numVertexies = 5;
        // pick two variables each time
        // and compute the cov value
        for (int p = 0; p < numVertexies; ++p)
        {
            for (int q = p; q < numVertexies; ++q)
            {
                float cov = find_covariance(inPointFieldVecEnsemble[p], inPointFieldVecEnsemble[q], inMeanArray[p], inMeanArray[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

        return;
    }

private:
    int m_NumSamples = 1000;
};

#endif //