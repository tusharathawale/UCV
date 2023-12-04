#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_EL_LESS_EIGENS_SANITY_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_EL_LESS_EIGENS_SANITY_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>

#include "./linalg/EasyLinalg/eigen.h"

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#else
// using the std library
#include <random>
#endif // VTKM_CUDA

class MVGaussianWithEnsemble2DTryELEntropyLessEigensSanity : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble2DTryELEntropyLessEigensSanity(double isovalue, int num_sample)
        : m_isovalue(isovalue), m_num_sample(num_sample){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell,
                                  FieldOutCell);
#ifdef DEBUG_WORKLET
    using ExecutionSignature = void(_2, _3, _4, _5, WorkIndex);
#else
    using ExecutionSignature = void(_2, _3, _4, _5);
#endif
    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble,
              typename OutCellFieldType1,
              typename OutCellFieldType2,
              typename OutCellFieldType3>

#ifdef DEBUG_WORKLET
    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy, vtkm::Id workIndex) const
#else
    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy) const
#endif
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numVertexies = inPointFieldVecEnsemble.GetNumberOfComponents();

        // TODO, using numVertexies to decide the length of mean and cov
        // and decide them at the runtime
        if (numVertexies != 4)
        {
            printf("the MVGaussianWithEnsemble2DTryLialg only support cell with 4 vertexies");
            return;
        }

        vtkm::Vec4f_64 meanArray;

        // get the type in the fieldVec
        // the VecType specifies the number of ensembles
        using VecType = decltype(inPointFieldVecEnsemble[0]);

        meanArray[0] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(0)]);
        meanArray[1] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(1)]);
        meanArray[2] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(2)]);
        meanArray[3] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(3)]);

        // debug input
#ifdef DEBUG_WORKLET
        if (workIndex == 50)
        {
            std::cout << "debug input" << std::endl;
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 15; j++)
                {
                    std::cout << inPointFieldVecEnsemble[i][j] << ",";
                }
                std::cout << std::endl;
            }
        }
#endif

        // if (fabs(meanArray[0]) < 0.000001 && fabs(meanArray[1]) < 0.000001 && fabs(meanArray[2]) < 0.000001 && fabs(meanArray[3]) < 0.000001)
        //{
        //     outCellFieldCProb = 0;
        //     return;
        // }

        // set the trim options to filter out values that does not contain the iso value
        // there is no cross prob for this values
        // find min and cell for all cell values
        double cellMin = vtkm::Infinity64();
        double cellMax = vtkm::NegativeInfinity64();
        for (int i = 0; i < 4; i++)
        {
            find_min_max<VecType>(inPointFieldVecEnsemble[i], cellMin, cellMax);
        }
#ifdef DEBUG_WORKLET
        if (workIndex == 50)
        {
            printf("---debug min %lf max %lf\n", cellMin, cellMax);
        }
#endif

        if (this->m_isovalue < cellMin || this->m_isovalue > cellMax)
        {
            outCellFieldCProb = 0;
            return;
        }

        // std::vector<double> cov_matrix;
        // for 4*4 matrix, there are 10 numbers at upper conner
        vtkm::Vec<vtkm::FloatDefault, 10> cov_matrix;
        vtkm::IdComponent index = 0;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                int updatep = updateIndex4(p);
                int updateq = updateIndex4(q);
                float cov = find_covariance<VecType>(inPointFieldVecEnsemble[updatep], inPointFieldVecEnsemble[updateq], meanArray[p], meanArray[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

        // generate sample
        EASYLINALG::Vec<double, 4> ucvmeanv;
        for (int i = 0; i < 4; i++)
        {
            ucvmeanv[i] = meanArray[i];
        }

        vtkm::IdComponent numSamples = m_num_sample;

        EASYLINALG::Matrix<double, 4, 4> ucvcov4by4;
        int covindex = 0;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                // use the elements at the top half
                // printf("%f ", cov_matrix[covindex]);
                ucvcov4by4[p][q] = cov_matrix[covindex];
                if (p != q)
                {
                    // assign value to another helf
                    ucvcov4by4[q][p] = ucvcov4by4[p][q];
                }
                covindex++;
            }
        }
#ifdef DEBUG_WORKLET
        if (workIndex == 50)
        {
            std::cout << "debug ucvcov4by4" << std::endl;
            ucvcov4by4.Show();
        }
#endif
        // EASYLINALG::Matrix<double, 4, 4> A = EASYLINALG::SymmEigenDecomposition(ucvcov4by4, 0.00001, 20);
        // Transform the iso value
        EASYLINALG::Vec<double, 4> transformIso(0);
        for (int i = 0; i < 4; i++)
        {
            transformIso[i] = this->m_isovalue - ucvmeanv[i];
        }
#ifdef DEBUG_WORKLET
        if (workIndex == 50)
        {
            std::cout << "debug transformIso " << std::endl;
            transformIso.Show();
        }
#endif
        // compute eigen values
        EASYLINALG::Vec<double, 4> eigenValues;
        EASYLINALG::SymmEigenValues(ucvcov4by4, this->m_tolerance, this->m_iterations, eigenValues);

        // the first index stores which it corresponds to
        // the second index stores the inner position of eigen value
        EASYLINALG::Vec<EASYLINALG::Vec<double, 4>, 4> eigenVectors;
        for (int i = 0; i < 4; i++)
        {
            eigenVectors[i] = EASYLINALG::ComputeEigenVectors(ucvcov4by4, eigenValues[i], this->m_iterations);
#ifdef DEBUG_WORKLET
            if (workIndex == 50)
            {
                std::cout << "debug eigen values " << eigenValues[i] << std::endl;
                eigenVectors[i].Show();
            }
#endif
        }
        // LIAG_FUNC_MACRO Vec<T, Size> ComputeEigenVectors(const Matrix<T, Size, Size> &A, const T &eigenValue, uint maxIter)

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng;
        thrust::random::normal_distribution<double> norm(0, 1);
#else
        std::mt19937 rng;
        rng.seed(std::mt19937::default_seed);
        std::normal_distribution<double> norm(0, 1);
#endif // VTKM_CUDA

        vtkm::Vec<vtkm::FloatDefault, 16> probHistogram;
        for (int i = 0; i < 16; i++)
        {
            probHistogram[i] = 0.0;
        }

        EASYLINALG::Vec<double, 4> sample_v;

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            // clear it sample results each time
            EASYLINALG::Vec<double, 4> sampleResults(0);
            // get sample vector
            // only use the largest eigen value currently
            // refer this for detailed ideas:
            // https://stephens999.github.io/fiveMinuteStats/mvnorm_eigen.html
            // only need to sample it one time
            for (int i = 0; i < 4; i++)
            {
                sample_v[i]=vtkm::Sqrt(eigenValues[i])*norm(rng);
                // std::normal_distribution<double> norm(0, vtkm::Sqrt(eigenValues[i]));
                // sample_v[i] = norm(rng);
            }

            // compute sampled results
            // for each sampled results
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    // be careful, each row in matrix is a eigen vector
                    sampleResults[i] += eigenVectors[j][i] * sample_v[j];
                    // sampleResults[i] += eigenVectors[j][i];
                }
            }
#ifdef DEBUG_WORKLET
            if (workIndex == 50 && n == 10)
            {
                std::cout << "debug sampleResults " << std::endl;
                sampleResults.Show();
            }
#endif
            // compute the specific position
            // map > or < to specific cases
            uint caseValue = 0;
            for (uint i = 0; i < 4; i++)
            {
                // setting associated position to 1 if iso is larger than specific cases
                if (transformIso[i] >= sampleResults[i])
                {
                    caseValue = (1 << i) | caseValue;
                }
                // the associated pos is 0 otherwise
            }
            probHistogram[caseValue] = probHistogram[caseValue] + 1.0;
        }

        // go through probHistogram and compute pro
        for (int i = 0; i < 16; i++)
        {
            probHistogram[i] = (probHistogram[i] / (1.0 * numSamples));
            // printf("debug caseValue %d probHistogram %f\n", i, probHistogram[i]);
        }

        // cross probability
        // outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);
        outCellFieldCProb = 1.0 - (probHistogram[0] + probHistogram[15]);

        vtkm::Id nonzeroCases = 0;
        vtkm::FloatDefault entropyValue = 0;
        vtkm::FloatDefault templog = 0;
        // compute number of nonzero cases
        // compute entropy
        for (int i = 0; i < 16; i++)
        {
            if (probHistogram[i] > 0.0001)
            {
                nonzeroCases++;
                templog = vtkm::Log2(probHistogram[i]);
                // if (i != 0 && i != totalNumCases - 1)
                //{
                //     totalnonzeroProb += probHistogram[i];
                // }
            }
            // do not update entropy if the pro is zero
            entropyValue = entropyValue + (-probHistogram[i]) * templog;
        }

        outCellFieldNumNonzeroProb = nonzeroCases;
        outCellFieldEntropy = entropyValue;
    }

    VTKM_EXEC int updateIndex4(int index) const
    {
        if (index == 0)
        {
            return 0;
        }
        else if (index == 1)
        {
            return 3;
        }
        else if (index == 2)
        {
            return 2;
        }
        else if (index == 3)
        {
            return 1;
        }

        printf("error, failed to compute updateIndex4\n");

        return 0;
    }

    template <typename VecType>
    VTKM_EXEC void find_min_max(const VecType &arr, vtkm::Float64 &min, vtkm::Float64 &max) const
    {
        vtkm::Id num = arr.GetNumberOfComponents();
        for (vtkm::Id i = 0; i < arr.GetNumberOfComponents(); i++)
        {
            min = vtkm::Min(min, arr[i]);
            max = vtkm::Max(max, arr[i]);
        }
        return;
    }

    template <typename VecType>
    VTKM_EXEC vtkm::Float64 find_mean(const VecType &arr) const
    {
        vtkm::Float64 sum = 0;
        vtkm::Id num = arr.GetNumberOfComponents();
        for (vtkm::Id i = 0; i < arr.GetNumberOfComponents(); i++)
        {
            sum = sum + arr[i];
        }
        vtkm::Float64 mean = (1.0 * sum) / (1.0 * num);
        return mean;
    }
    template <typename VecType>
    VTKM_EXEC double find_covariance(const VecType &arr1, const VecType &arr2,
                                     double &mean1, double &mean2) const
    {
        if (arr1.GetNumberOfComponents() != arr2.GetNumberOfComponents())
        {
            // cuda does not support exception
            printf("error, failed to compute find_covariance, the array size should be equal with each other\n");
            return 0;
        }
        vtkm::Id arraySize = arr1.GetNumberOfComponents();
        double sum = 0;
        for (int i = 0; i < arraySize; i++)
            sum = sum + (arr1[i] - mean1) * (arr2[i] - mean2);
        return (double)sum / (double)(arraySize - 1);
    }

private:
    double m_isovalue;
    int m_num_sample = 1000;
    int m_iterations = 200;
    double m_tolerance = 0.00001;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h
