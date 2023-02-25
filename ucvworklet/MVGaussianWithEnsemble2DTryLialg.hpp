#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>

// #include "./linalg/ucv_matrix.h"
#include "./linalg/ucv_matrix_static_4by4.h"

// use this as the results checking on cpu
// #include "./eigenmvn.h"
// this worklet is for the input data that put the different data in a separate array
// for the wind data here https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200
// there are 15 numbers (ensemble extraction) each data is put in a different file

class MVGaussianWithEnsemble2DTryLialg : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble2DTryLialg(double isovalue, int num_sample)
        : m_isovalue(isovalue), m_num_sample(num_sample){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, WorkIndex);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble, typename OutCellFieldType>

    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutCellFieldType &outCellFieldCProb, vtkm::Id workIndex) const
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

        // set the trim options to filter the 0 values
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

        UCVMATH::vec_t ucvmeanv;

        for (int i = 0; i < 4; i++)
        {
            ucvmeanv.v[i] = meanArray[i];
        }

        vtkm::IdComponent numSamples = m_num_sample;
        vtkm::Id numCrossings = 0;
        // this can be adapted to 3d case

        UCVMATH::mat_t ucvcov4by4;
        int covindex = 0;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                // use the elements at the top half
                // printf("%f ", cov_matrix[covindex]);
                ucvcov4by4.v[p][q] = cov_matrix[covindex];
                if (p != q)
                {
                    // assign value to another helf
                    ucvcov4by4.v[q][p] = ucvcov4by4.v[p][q];
                }
                covindex++;
            }
        }

        // if (workIndex ==15822)
        //{
        //     matrix_show(&ucvcov4by4);
        // }

        double result[4];
        eigen_solve_eigenvalues(&ucvcov4by4, 0.000001, 50, result);

        UCVMATH::mat_t A = UCVMATH::eigen_vector_decomposition(&ucvcov4by4);

        // if (workIndex ==15822)
        //{
        //     printf("eigen values\n");
        //     printf("%8.7f %8.7f %8.7f %8.7f\n",result[0],result[1],result[2],result[3]);
        //     printf("eigen dec ucv\n");
        //     matrix_show(&A);
        // }

        /*
        using the eigen to do the same operation to check the results

        Eigen::Vector4d meanVector(4);
        meanVector << meanArray[0], meanArray[1], meanArray[2], meanArray[3];
        Eigen::Matrix4d cov4by4(4, 4);
        cov4by4 << cov_matrix[0], cov_matrix[1], cov_matrix[2], cov_matrix[3],
            cov_matrix[1], cov_matrix[4], cov_matrix[5], cov_matrix[6],
            cov_matrix[2], cov_matrix[5], cov_matrix[7], cov_matrix[8],
            cov_matrix[3], cov_matrix[6], cov_matrix[8], cov_matrix[9];
        Eigen::EigenMultivariateNormal<vtkm::Float64> normX_solver(meanVector, cov4by4);
        if (workIndex ==15822)
        {
            std::cout << normX_solver._eigenSolver.eigenvalues() << std::endl;
            printf("eigen dec, eigen lib\n");
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    printf(" %f", normX_solver._transform(i, j));
                }
                printf("\n");
            }
        }
        */

        UCVMATH::vec_t sample_v;
        UCVMATH::vec_t AUM;

#if defined(VTKM_CUDA) || defined(VTKM_KOKKOS_HIP)
        thrust::minstd_rand rng;
        thrust::random::normal_distribution<double> norm;
#else
        std::mt19937 rng;
        rng.seed(std::mt19937::default_seed);
        std::normal_distribution<double> norm;
#endif // VTKM_CUDA

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            // get sample vector
            for (int i = 0; i < 4; i++)
            {
                // using other sample mechanism such as thrust as needed
                sample_v.v[i] = norm(rng);
            }

            AUM = UCVMATH::matrix_mul_vec_add_vec(&A, &sample_v, &ucvmeanv);

            // if (n==0 && workIndex == 15822)
            //{
            //     printf("ucv AUM\n");
            //     vec_show(&AUM);
            // }

            if ((m_isovalue <= AUM.v[0]) && (m_isovalue <= AUM.v[1]) && (m_isovalue <= AUM.v[2]) && (m_isovalue <= AUM.v[3]))
            {
                numCrossings = numCrossings + 0;
            }
            else if ((m_isovalue >= AUM.v[0]) && (m_isovalue >= AUM.v[1]) && (m_isovalue >= AUM.v[2]) && (m_isovalue >= AUM.v[3]))
            {
                numCrossings = numCrossings + 0;
            }
            else
            {
                numCrossings = numCrossings + 1;
            }
        }

        // cross probability
        // std::cout << "ucv numCrossings " << numCrossings << std::endl;
        outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);
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
