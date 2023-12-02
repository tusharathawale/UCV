#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_EL_EIGEN_DECOMP_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_EL_EIGEN_DECOMP_h

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

class MVGaussianWithEnsemble2DTryELEigenDecomp : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble2DTryELEigenDecomp(int num_sample)
        : m_num_sample(num_sample){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble,
              typename OutMatrixType,
              typename OutVecType>

    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutMatrixType &outEigenDecompMatrix,
        OutVecType &outMeanArray) const
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

        // get the type in the fieldVec
        // the VecType specifies the number of ensembles
        using VecType = decltype(inPointFieldVecEnsemble[0]);

        outMeanArray[0] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(0)]);
        outMeanArray[1] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(1)]);
        outMeanArray[2] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(2)]);
        outMeanArray[3] = find_mean<VecType>(inPointFieldVecEnsemble[updateIndex4(3)]);

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
                float cov = find_covariance<VecType>(inPointFieldVecEnsemble[updatep], inPointFieldVecEnsemble[updateq], outMeanArray[p], outMeanArray[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

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

        EASYLINALG::Matrix<double, 4, 4> A = EASYLINALG::SymmEigenDecomposition(ucvcov4by4, this->m_tolerance, this->m_iterations);

        for (int j = 0; j < 4; j++)
        {
            for (int i = 0; i < 4; i++)
            {
                outEigenDecompMatrix[j][i] = A[j][i];
            }
        }
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
    int m_num_sample = 1000;
    int m_iterations = 200;
    double m_tolerance = 0.00001;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h
