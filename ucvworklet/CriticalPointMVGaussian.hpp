#ifndef UCV_CRITICAL_POINT_MV_GAUSSIAN_h
#define UCV_CRITICAL_POINT_MV_GAUSSIAN_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>
#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
#include <thrust/device_vector.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#else
#include <random>
#endif

#include "./linalg/EasyLinalg/eigen.h"

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

    using ExecutionSignature = void(_2, _3, Boundary, WorkIndex);

    template <typename InPointField, typename OutPointField>
    VTKM_EXEC void operator()(const InPointField &inPointFieldVecEnsemble,
                              OutPointField &minProb,
                              const vtkm::exec::BoundaryState &boundary,
                              vtkm::Id workIndex) const
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

        using VecType = decltype(inPointFieldVecEnsemble.Get(0, 0, 0));
        // using VecType2 = decltype(inPointFieldVecEnsemble);

        // vtkm::FloatDefault m1 = minValue.Get(0, 0, 0);
        // VecType ensArray0 = getEns<InPointField, VecType>(inPointFieldVecEnsemble, 0);
        // if (workIndex == 0)
        // {
        //     // maybe double look at how
        //     // printSummary_ArrayHandle get the actual type
        //     // printf("--%s\n", typeid(VecType2).name());
        //     printf("--get number of component %d\n", ensArray0.GetNumberOfComponents());
        // }

        // filter out the element in the boundry
        // if the element is at the boundry places its min prob is 0
        if ((maxIndices[0] - minIndices[0] < 2) || (maxIndices[1] - minIndices[1] < 2))
        {
            // if x and y is at the boundry, do not consider it
            minProb = 0;
            return;
        }

        // compute mean
        EASYLINALG::Vec<double, 5> ucvmeanv;

        for (vtkm::Id i = 0; i < 5; i++)
        {
            ucvmeanv[i] = find_mean<VecType>(getEns<InPointField, VecType>(inPointFieldVecEnsemble, i));
        }

        // filter out empty part
        if (fabs(ucvmeanv[0]) < 0.0000001 && fabs(ucvmeanv[1]) < 0.0000001 && fabs(ucvmeanv[2]) < 0.0000001 && fabs(ucvmeanv[3]) < 0.0000001 && fabs(ucvmeanv[4]) < 0.0000001)
        {
            minProb = 0;
            return;
        }

        vtkm::Vec<vtkm::Float64, 15> cov_matrix;
        vtkm::IdComponent index = 0;
        for (int p = 0; p < 5; ++p)
        {
            for (int q = p; q < 5; ++q)
            {
                float cov = find_covariance<VecType>(getEns<InPointField, VecType>(inPointFieldVecEnsemble, p), getEns<InPointField, VecType>(inPointFieldVecEnsemble, q), ucvmeanv[p], ucvmeanv[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

        // put cov matrix value into the EASYLINALG
        EASYLINALG::Matrix<double, 5, 5> ucvcov5by5;
        int covindex = 0;
        for (int p = 0; p < 5; ++p)
        {
            for (int q = p; q < 5; ++q)
            {
                // use the elements at the top half
                // printf("%f ", cov_matrix[covindex]);
                ucvcov5by5[p][q] = cov_matrix[covindex];

                if (p != q)
                {
                    // assign value to another helf
                    ucvcov5by5[q][p] = ucvcov5by5[p][q];
                }
                covindex++;
            }
        }

        EASYLINALG::Matrix<double, 5, 5> A = EASYLINALG::SymmEigenDecomposition(ucvcov5by5, this->m_tolerance, this->m_iterations);
        EASYLINALG::Vec<double, 5> sample_v;
        EASYLINALG::Vec<double, 5> AUM;

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng;
        thrust::random::normal_distribution<double> norm;
#else
        std::mt19937 rng;
        rng.seed(std::mt19937::default_seed);
        std::normal_distribution<double> norm;
#endif // VTKM_CUDA
        double localminCases = 0;
        for (vtkm::Id n = 0; n < this->m_NumSamples; ++n)
        {
            // get sample vector
            for (int i = 0; i < 5; i++)
            {
                sample_v[i] = norm(rng);
            }

            // Ax+b operation
            AUM = EASYLINALG::DGEMV(1.0, A, sample_v, 1.0, ucvmeanv);

            // TODO, compute minprob according to AUM
            // if x1 < x2,x1<x3 ... x1<x5
            if (AUM[0] < AUM[1] && AUM[0] < AUM[2] && AUM[0] < AUM[3] && AUM[0] < AUM[4])
            {
                localminCases = localminCases + 1;
            }
        }
        minProb = localminCases / this->m_NumSamples;

        return;
    }

    template <typename VecType>
    VTKM_EXEC vtkm::Float64 find_mean(const VecType &arr) const
    {
        vtkm::Float64 sum = 0;
        vtkm::Id num = arr.GetNumberOfComponents();
        for (vtkm::Id i = 0; i < arr.GetNumberOfComponents(); i++)
        {
            vtkm::Float64 temp = arr[i];
            sum = sum + temp;
        }
        vtkm::Float64 mean = (1.0 * sum) / (1.0 * num);
        return mean;
    }

    template <typename VecType>
    VTKM_EXEC vtkm::Float64 find_covariance(const VecType &arr1, const VecType &arr2,
                                            vtkm::Float64 &mean1, vtkm::Float64 &mean2) const
    {
        if (arr1.GetNumberOfComponents() != arr2.GetNumberOfComponents())
        {
            // cuda does not support exception
            printf("error, failed to compute find_covariance, the array size should be equal with each other\n");
            return 0;
        }
        vtkm::Id arraySize = arr1.GetNumberOfComponents();
        vtkm::Float64 sum = 0;
        for (int i = 0; i < arraySize; i++)
        {
            vtkm::Float64 temp1 = arr1[i];
            vtkm::Float64 temp2 = arr2[i];
            sum = sum + (temp1 - mean1) * (temp2 - mean2);
        }
        return (vtkm::Float64)sum / (vtkm::Float64)(arraySize - 1);
    }

    template <typename InPointField, typename EnsArrayType>
    VTKM_EXEC inline EnsArrayType getEns(const InPointField &inPointFieldVecEnsemble, vtkm::Id index) const
    {
        if (index == 0)
        {
            return inPointFieldVecEnsemble.Get(0, 0, 0);
        }
        else if (index == 1)
        {
            return inPointFieldVecEnsemble.Get(0, 1, 0);
        }
        else if (index == 2)
        {
            return inPointFieldVecEnsemble.Get(0, -1, 0);
        }
        else if (index == 3)
        {
            return inPointFieldVecEnsemble.Get(1, 0, 0);
        }
        else if (index == 4)
        {
            return inPointFieldVecEnsemble.Get(-1, 0, 0);
        }

        printf("Error, fail to getEns by index %d\n", index);
        return inPointFieldVecEnsemble.Get(0, 0, 0);
    }

private:
    int m_NumSamples = 1000;
    int m_iterations = 200;
    double m_tolerance = 0.00001;
};

#endif //