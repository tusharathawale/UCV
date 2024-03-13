#ifndef UCV_CRITICAL_POINT_MONTE_CARLO_h
#define UCV_CRITICAL_POINT_MONTE_CARLO_h

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

struct CriticalPointMonteCarloWorklet : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    CriticalPointMonteCarloWorklet(vtkm::Id samples) : m_NumSamples(samples){};

    using ControlSignature = void(CellSetIn, FieldInNeighborhood, FieldInNeighborhood, FieldOut);

    using ExecutionSignature = void(_2, _3, _4, Boundary, WorkIndex);

    template <typename InPointField, typename OutPointField>
    VTKM_EXEC void operator()(const InPointField &minValue,
                              const InPointField &maxValue,
                              OutPointField &minProb,
                              const vtkm::exec::BoundaryState &boundary,
                              vtkm::Id WorkIndex) const
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

        vtkm::FloatDefault a1 = minValue.Get(0, 0, 0);
        vtkm::FloatDefault b1 = maxValue.Get(0, 0, 0);

        vtkm::FloatDefault a2 = minValue.Get(0, 1, 0);
        vtkm::FloatDefault b2 = maxValue.Get(0, 1, 0);

        vtkm::FloatDefault a3 = minValue.Get(0, -1, 0);
        vtkm::FloatDefault b3 = maxValue.Get(0, -1, 0);

        vtkm::FloatDefault a4 = minValue.Get(1, 0, 0);
        vtkm::FloatDefault b4 = maxValue.Get(1, 0, 0);

        vtkm::FloatDefault a5 = minValue.Get(-1, 0, 0);
        vtkm::FloatDefault b5 = maxValue.Get(-1, 0, 0);

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng;
        thrust::uniform_real_distribution<vtkm::FloatDefault> GenerateV1(a1, b1);
        thrust::uniform_real_distribution<vtkm::FloatDefault> GenerateV2(a2, b2);
        thrust::uniform_real_distribution<vtkm::FloatDefault> GenerateV3(a3, b3);
        thrust::uniform_real_distribution<vtkm::FloatDefault> GenerateV4(a4, b4);
        thrust::uniform_real_distribution<vtkm::FloatDefault> GenerateV5(a5, b5);
        vtkm::Id NumMinCase = 0;
        for (vtkm::Id i = 0; i < this->m_NumSamples; i++)
        {
            vtkm::FloatDefault V1 = GenerateV1(rng);
            vtkm::FloatDefault V2 = GenerateV2(rng);
            vtkm::FloatDefault V3 = GenerateV3(rng);
            vtkm::FloatDefault V4 = GenerateV4(rng);
            vtkm::FloatDefault V5 = GenerateV5(rng);

            if (V1 < V2 && V1 < V3 && V1 < V4 && V1 < V5)
            {
                NumMinCase++;
            }
        }
#else
        std::random_device rd;
        std::mt19937 gen(rd());

        // sample points for each range
        //  Generate samples from data rectangle
        std::uniform_real_distribution<vtkm::FloatDefault> GenerateV1(a1, b1);
        std::uniform_real_distribution<vtkm::FloatDefault> GenerateV2(a2, b2);
        std::uniform_real_distribution<vtkm::FloatDefault> GenerateV3(a3, b3);
        std::uniform_real_distribution<vtkm::FloatDefault> GenerateV4(a4, b4);
        std::uniform_real_distribution<vtkm::FloatDefault> GenerateV5(a5, b5);

        vtkm::Id NumMinCase = 0;
        for (vtkm::Id i = 0; i < this->m_NumSamples; i++)
        {
            vtkm::FloatDefault V1 = GenerateV1(gen);
            vtkm::FloatDefault V2 = GenerateV2(gen);
            vtkm::FloatDefault V3 = GenerateV3(gen);
            vtkm::FloatDefault V4 = GenerateV4(gen);
            vtkm::FloatDefault V5 = GenerateV5(gen);

            if (V1 < V2 && V1 < V3 && V1 < V4 && V1 < V5)
            {
                NumMinCase++;
            }
        }

#endif

        // compute the min prob
        minProb = (NumMinCase * 1.0) / (1.0 * m_NumSamples);
        return;
    }

private:
    int m_NumSamples = 1000;
};

#endif // UCV_CRITICAL_POINT_h