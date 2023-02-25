#ifndef UCV_MULTIVARIANT_GAUSSIAN3D_h
#define UCV_MULTIVARIANT_GAUSSIAN3D_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>
//#include <Eigen/Dense>
#include "./linalg/ucv_matrix_static_8by8.h"

class MVGaussianWithEnsemble3DTryLialg : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble3DTryLialg(double isovalue, int numSamples)
        : m_isovalue(isovalue), m_numSamples(numSamples){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell,
                                  FieldOutCell,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4, _5, _6);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble,
              typename InPointFieldVecMean,
              typename OutCellFieldType1,
              typename OutCellFieldType2,
              typename OutCellFieldType3>
    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        const InPointFieldVecMean &inMeanArray,
        OutCellFieldType1 &outCellFieldCProb,
        OutCellFieldType2 &outCellFieldNumNonzeroProb,
        OutCellFieldType3 &outCellFieldEntropy) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numVertexies = inPointFieldVecEnsemble.GetNumberOfComponents();
        const uint8_t numVertex3d = 8;
        if (numVertexies != numVertex3d)
        {
            // throw std::runtime_error("MVGaussianWithEnsemble3DTryLialg expects 8 vertecies");
            printf("MVGaussianWithEnsemble3DTryLialg expects 8 vertecies\n");
            return;
        }

        if (inMeanArray.GetNumberOfComponents() != numVertex3d)
        {
            // throw std::runtime_error("inMeanArray in MVGaussianWithEnsemble3DTryLialg expects 8 vertecies");
            printf("inMeanArray in MVGaussianWithEnsemble3DTryLialg expects 8 vertecies\n");
            return;
        }

        if (inPointFieldVecEnsemble[0].GetNumberOfComponents() != numVertex3d * numVertex3d)
        {
            // throw std::runtime_error("only support ensemble size 64 for blockSize equals to 4");
            printf("only support ensemble size 64 for blockSize equals to 4\n");
            return;
        }

        // std::vector<double> cov_matrix;
        vtkm::Vec<vtkm::FloatDefault, 36> cov_matrix;
        vtkm::IdComponent index = 0;
        for (int p = 0; p < numVertexies; ++p)
        {
            for (int q = p; q < numVertexies; ++q)
            {
                float cov = find_covariance(inPointFieldVecEnsemble[p], inPointFieldVecEnsemble[q], inMeanArray[p], inMeanArray[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

        // generate sample
        UCVMATH::vec_t ucvmeanv;

        for (int i = 0; i < numVertex3d; i++)
        {
            ucvmeanv.v[i] = inMeanArray[i];
        }

        // generate mean and cov matrix
        UCVMATH::mat_t ucvcov8by8;
        int covindex = 0;
        for (int p = 0; p < numVertex3d; ++p)
        {
            for (int q = p; q < numVertex3d; ++q)
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

        vtkm::IdComponent numSamples = this->m_numSamples;

        UCVMATH::mat_t A = UCVMATH::eigen_vector_decomposition(&ucvcov8by8);

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

        vtkm::Vec<vtkm::FloatDefault, 256> probHistogram;

        // init to 0
        for (int i = 0; i < 256; i++)
        {
            probHistogram[i] = 0.0;
        }

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            // std::cout << R.coeff(n, 0) << " " << R.coeff(n, 1) << " " << R.coeff(n, 2) << " " << R.coeff(n, 3) << std::endl;
            // get sample vector
            for (int i = 0; i < numVertex3d; i++)
            {
                // using other sample mechanism such as thrust as needed
                sample_v.v[i] = norm(rng);
            }

            AUM = UCVMATH::matrix_mul_vec_add_vec(&A, &sample_v, &ucvmeanv);

            /*
            if ((m_isovalue <= AUM.v[0]) && (m_isovalue <= AUM.v[1]) && (m_isovalue <= AUM.v[2]) && (m_isovalue <= AUM.v[3]) && (m_isovalue <= AUM.v[4]) && (m_isovalue <= AUM.v[5]) && (m_isovalue <= AUM.v[6]) && (m_isovalue <= AUM.v[7]))
            {
                numCrossings = numCrossings + 0;
            }
            else if ((m_isovalue >= AUM.v[0]) && (m_isovalue >= AUM.v[1]) && (m_isovalue >= AUM.v[2]) && (m_isovalue >= AUM.v[3]) && (m_isovalue >= AUM.v[4]) && (m_isovalue >= AUM.v[5]) && (m_isovalue >= AUM.v[6]) && (m_isovalue >= AUM.v[7]))
            {
                numCrossings = numCrossings + 0;
            }
            else
            {
                //there are 254 total cases here 256 - all 0 cases - all 1 cases
                numCrossings = numCrossings + 1;
            }
            */

            // go through 8 cases
            uint caseValue = 0;
            for (uint i = 0; i < 8; i++)
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
        for (int i = 0; i < 256; i++)
        {
            probHistogram[i] = (probHistogram[i] / (1.0 * numSamples));
            //printf("debug caseValue %d probHistogram %f\n", i, probHistogram[i]);
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

    // how to get this vtkm::Vec<double, 64> in an more efficient way
    vtkm::FloatDefault find_mean(const vtkm::Vec<double, 64> &arr) const
    {
        vtkm::FloatDefault sum = 0;
        vtkm::Id num = arr.GetNumberOfComponents();
        for (vtkm::Id i = 0; i < arr.GetNumberOfComponents(); i++)
        {
            sum = sum + arr[i];
        }
        vtkm::FloatDefault mean = (vtkm::FloatDefault)sum / (vtkm::FloatDefault)(num);
        return mean;
    }

    VTKM_EXEC inline vtkm::FloatDefault find_covariance(const vtkm::Vec<vtkm::FloatDefault, 64> &arr1, const vtkm::Vec<vtkm::FloatDefault, 64> &arr2,
                                                        const vtkm::FloatDefault &mean1, const vtkm::FloatDefault &mean2) const
    {
        if (arr1.GetNumberOfComponents() != arr2.GetNumberOfComponents())
        {
            printf("error, failed to compute find_covariance, the array size should be equal with each other\n");
            return 0;
        }
        vtkm::Id arraySize = arr1.GetNumberOfComponents();
        vtkm::FloatDefault sum = 0;
        for (int i = 0; i < arraySize; i++)
            sum = sum + (arr1[i] - mean1) * (arr2[i] - mean2);
        return sum / (vtkm::FloatDefault)(arraySize - 1);
    }

private:
    double m_isovalue;
    int m_numSamples;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN3D_h