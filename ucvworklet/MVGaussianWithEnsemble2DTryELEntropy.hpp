#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_EL_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_EL_h

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

class MVGaussianWithEnsemble2DTryELEntropy : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble2DTryELEntropy(double isovalue, int num_sample)
        : m_isovalue(isovalue), m_num_sample(num_sample){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4, _5, WorkIndex);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble,
              typename OutCellFieldType1,
              typename OutCellFieldType2,
              typename OutCellFieldType3>

    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy,
        vtkm::Id workIndex) const
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

        // set the trim options to filter out the 0 values
        if (fabs(meanArray[0]) < 0.000001 && fabs(meanArray[1]) < 0.000001 && fabs(meanArray[2]) < 0.000001 && fabs(meanArray[3]) < 0.000001)
        {
            outCellFieldCProb = 0;
            return;
        }

        // if (workIndex == 0)
        //{
        //     std::cout << meanArray[0] << " " << meanArray[1] << " " << meanArray[2] << " " << meanArray[3] << std::endl;
        // }

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

        // UCVMATH::vec_t ucvmeanv;
        // gsl_vector *ucvmeanv = UCVMATH_CSTM_GSL::cstm_gsl_vector_alloc(4);
        EASYLINALG::Vec<double, 4> ucvmeanv;
        for (int i = 0; i < 4; i++)
        {
            ucvmeanv[i] = meanArray[i];
        }

        vtkm::IdComponent numSamples = m_num_sample;
        // vtkm::Id numCrossings = 0;
        // this can be adapted to 3d case

        // UCVMATH::mat_t ucvcov4by4_original;
        // gsl_matrix *ucvcov4by4 = gsl_matrix_alloc(4, 4);
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

        // if (workIndex ==15822)
        //{
        //     matrix_show(&ucvcov4by4);
        // }

        // UCVMATH::mat_t AOriginal = UCVMATH::eigen_vector_decomposition(&ucvcov4by4_original);
        // gsl_matrix *A = UCVMATH_CSTM_GSL::gsl_eigen_vector_decomposition(ucvcov4by4);
        EASYLINALG::Matrix<double, 4, 4> A = EASYLINALG::SymmEigenDecomposition(ucvcov4by4, 0.00001, 20);

        // some values are filtered out since it can be in the empty region
        // with 0 values there
        
        if (workIndex == 9896)
        {
            printf("index is %d\n",workIndex);
            printf("matrix ucvcov4by4\n");
            ucvcov4by4.Show();
            printf("matrix A\n");
            A.Show();
        }
        

        // UCVMATH::vec_t sample_v;
        // UCVMATH::vec_t AUM;
        // gsl_vector *sample_v = UCVMATH_CSTM_GSL::cstm_gsl_vector_alloc(4);
        // gsl_vector *AUM = UCVMATH_CSTM_GSL::cstm_gsl_vector_alloc(4);
        EASYLINALG::Vec<double, 4> sample_v;
        EASYLINALG::Vec<double, 4> AUM;

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng;
        thrust::random::normal_distribution<double> norm;
#else
        std::mt19937 rng;
        rng.seed(std::mt19937::default_seed);
        std::normal_distribution<double> norm;
#endif // VTKM_CUDA

        vtkm::Vec<vtkm::FloatDefault, 16> probHistogram;
        for (int i = 0; i < 16; i++)
        {
            probHistogram[i] = 0.0;
        }

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            // get sample vector
            for (int i = 0; i < 4; i++)
            {
                // using other sample mechanism such as thrust as needed
                // sample_v.v[i] = norm(rng);
                // gsl_vector_set(sample_v,i,norm(rng));
                sample_v[i] = norm(rng);
            }

            // Ax+b operation
            AUM = EASYLINALG::DGEMV(1.0, A, sample_v, 1.0, ucvmeanv);

            // compute the specific position
            // map > or < to specific cases
            uint caseValue = 0;
            for (uint i = 0; i < 4; i++)
            {
                // setting associated position to 1 if iso larger then specific cases
                if (m_isovalue >= AUM[i])
                {
                    caseValue = (1 << i) | caseValue;
                }
            }

            // the associated pos is 0 otherwise
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
            return 1;
        }
        else if (index == 3)
        {
            return 2;
        }

        printf("error, failed to compute updateIndex4\n");

        return 0;
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
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h
