#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_EL_SAMPLING_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_EL_SAMPLING_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/cont/ArrayHandleRandomStandardNormal.h>

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#else
// using the std library
#include <random>
#endif // VTKM_CUDA

class MVGaussianWithEnsemble2DSampling : public vtkm::worklet::WorkletMapField
{
public:
    MVGaussianWithEnsemble2DSampling(){};

    using ControlSignature = void(FieldInOut);

    using ExecutionSignature = void(_1,WorkIndex);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename ValueType>
    VTKM_EXEC void operator()(
        ValueType &sampleVec,
        vtkm::Id workIndex) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numValues = sampleVec.GetNumberOfComponents();
        #if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng(workIndex);
        thrust::random::normal_distribution<double> norm;
#else
        std::mt19937 rng;
        //need to set offset here, otherwise, the seeds are same for all thread
        rng.seed(std::mt19937::default_seed+workIndex);
        std::normal_distribution<double> norm;
#endif // VTKM_CUDA
        for (int i = 0; i < numValues; i++)
        {
            sampleVec[i] = norm(rng);
            // for debug
            // sampleVec[i] = 0.1;
        }
    }
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h
