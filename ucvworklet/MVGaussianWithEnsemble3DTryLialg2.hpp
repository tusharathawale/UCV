#ifndef UCV_MULTIVARIANT_GAUSSIAN3D2_h
#define UCV_MULTIVARIANT_GAUSSIAN3D2_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>

// #include "./linalg/ucv_matrix.h"
#include "./linalg/ucv_matrix_static_8by8.h"

class MVGaussianWithEnsemble3DTryLialg2 : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble3DTryLialg2(double isovalue, int num_sample)
        : m_isovalue(isovalue), m_num_sample(num_sample){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4, _5);

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
        OutCellFieldType3 &outCellFieldEntropy) const
    {
        // TODO try cuda function
        // get thread block number
        // grid size etc 

        // how to process the case where there are multiple variables
        vtkm::IdComponent numVertexies = inPointFieldVecEnsemble.GetNumberOfComponents();

        // TODO, using numVertexies to decide the length of mean and cov
        // and decide them at the runtime
        if (numVertexies != 8)
        {
            printf("the MVGaussianWithEnsemble3DTryLialg2 only support cell with 4 vertexies");
            return;
        }

        // get the type in the fieldVec
        // the VecType specifies the number of ensembles
        using VecType = decltype(inPointFieldVecEnsemble[0]);
        UCVMATH::vec_t ucvmeanv;

        for (int i = 0; i < numVertexies; i++)
        {
            ucvmeanv.v[i] = find_mean<VecType>(inPointFieldVecEnsemble[i]);
        }

        // set the trim options to filter the 0 values
        if (fabs(ucvmeanv.v[0]) < 0.000001 && fabs(ucvmeanv.v[1]) < 0.000001 && fabs(ucvmeanv.v[2]) < 0.000001 && fabs(ucvmeanv.v[3]) < 0.000001 && fabs(ucvmeanv.v[4]) < 0.000001 && fabs(ucvmeanv.v[5]) < 0.000001 && fabs(ucvmeanv.v[6]) < 0.000001 && fabs(ucvmeanv.v[7]) < 0.000001)
        {
            outCellFieldCProb = 0;
            return;
        }

        // if (workIndex == 0)
        //{
        //     std::cout << meanArray[0] << " " << meanArray[1] << " " << meanArray[2] << " " << meanArray[3] << std::endl;
        // }

        // std::vector<double> cov_matrix;
        // for 8*8 matrix, there are 36 numbers at upper conner
        vtkm::Vec<vtkm::FloatDefault, 36> cov_matrix;
        vtkm::IdComponent index = 0;
        for (int p = 0; p < numVertexies; ++p)
        {
            for (int q = p; q < numVertexies; ++q)
            {
                float cov = find_covariance<VecType>(inPointFieldVecEnsemble[p], inPointFieldVecEnsemble[q], ucvmeanv.v[p], ucvmeanv.v[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

        // generate sample

        vtkm::IdComponent numSamples = m_num_sample;
        // vtkm::Id numCrossings = 0;
        // this can be adapted to 3d case

        UCVMATH::mat_t ucvcov8by8;
        int covindex = 0;
        for (int p = 0; p < numVertexies; ++p)
        {
            for (int q = p; q < numVertexies; ++q)
            {
                // use the elements at the top half
                // printf("%f ", cov_matrix[covindex]);
                ucvcov8by8.v[p][q] = cov_matrix[covindex];
                if (p != q)
                {
                    // assign value to another helf
                    ucvcov8by8.v[q][p] = ucvcov8by8.v[p][q];
                }
                covindex++;
            }
        }

        // if (workIndex ==15822)
        //{
        //     matrix_show(&ucvcov4by4);
        // }

        //double result[8];
        //eigen_solve_eigenvalues(&ucvcov8by8, 0.000001, 50, result);

        UCVMATH::mat_t A = UCVMATH::eigen_vector_decomposition(&ucvcov8by8);

        UCVMATH::vec_t sample_v;
        UCVMATH::vec_t AUM;

#ifdef VTKM_CUDA
        thrust::minstd_rand rng;
        thrust::random::normal_distribution<double> norm;
#else
        std::mt19937 rng;
        rng.seed(std::mt19937::default_seed);
        std::normal_distribution<double> norm;
#endif // VTKM_CUDA

        vtkm::Vec<vtkm::FloatDefault, 256> probHistogram;
        for (int i = 0; i < 256; i++)
        {
            probHistogram[i] = 0.0;
        }

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            // get sample vector
            for (int i = 0; i < numVertexies; i++)
            {
                // using other sample mechanism such as thrust as needed
                sample_v.v[i] = norm(rng);
            }

            AUM = UCVMATH::matrix_mul_vec_add_vec(&A, &sample_v, &ucvmeanv);

            // compute the specific position
            // map > or < to specific cases
            uint caseValue = 0;
            for (uint i = 0; i < numVertexies; i++)
            {
                // setting associated position to 1 if iso larger then specific cases
                if (m_isovalue >= AUM.v[i])
                {
                    caseValue = (1 << i) | caseValue;
                }
            }

            // the associated pos is 0 otherwise
            probHistogram[caseValue] = probHistogram[caseValue] + 1.0;
        }

        // go through probHistogram and compute pro
        for (int i = 0; i < numVertexies*numVertexies; i++)
        {
            probHistogram[i] = (probHistogram[i] / (1.0 * numSamples));
            // printf("debug caseValue %d probHistogram %f\n", i, probHistogram[i]);
        }

        // cross probability
        // outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);
        outCellFieldCProb = 1.0 - (probHistogram[0] + probHistogram[255]);

        vtkm::Id nonzeroCases = 0;
        vtkm::FloatDefault entropyValue = 0;
        vtkm::FloatDefault templog = 0;
        // compute number of nonzero cases
        // compute entropy
        for (int i = 0; i < 256; i++)
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
