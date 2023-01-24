#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>

#include "./linalg/ucv_matrix.h"

// use this as the results checking on cpu
// #include "./eigenmvn.h"
// this worklet is for the input data that put the different data in a separate array
// for the wind data here https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200
// there are 15 numbers (ensemble extraction) each data is put in a different file

class MVGaussianWithEnsemble2DTryLialg : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble2DTryLialg(double isovalue)
        : m_isovalue(isovalue){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble, typename OutCellFieldType>

    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        OutCellFieldType &outCellFieldCProb) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numVertexies = inPointFieldVecEnsemble.GetNumberOfComponents();

        // TODO, extracting data from 4 vertexies and compute the uncertainty things
        // this is supposed to be 4
        // std::cout << "size of numVertexies " << numVertexies << std::endl;
        // this is supposed to be 15
        // std::cout << "ensemble number " << inPointFieldVecEnsemble[0].GetNumberOfComponents() << std::endl;

        // this InPointFieldVecEnsemble here is supposed to be the Vec15
        // TODO, compute the mean, cov and cross probability

        // vector is not good in cuda
        // std::vector<vtkm::Float64> meanArray(4, 0);
        vtkm::Vec<vtkm::FloatDefault, 4> meanArray;
        // derive type
        meanArray[0] = find_mean(inPointFieldVecEnsemble[0]);
        meanArray[1] = find_mean(inPointFieldVecEnsemble[1]);
        meanArray[2] = find_mean(inPointFieldVecEnsemble[2]);
        meanArray[3] = find_mean(inPointFieldVecEnsemble[3]);

        // std::vector<double> cov_matrix;
        // for 4*4 matrix, there are 10 numbers at upper conner
        vtkm::Vec<vtkm::FloatDefault, 10> cov_matrix;
        vtkm::IdComponent index = 0;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                float cov = find_covariance(inPointFieldVecEnsemble[p], inPointFieldVecEnsemble[q], meanArray[p], meanArray[q]);
                cov_matrix[index] = cov;
                index++;
            }
        }

        // generate sample

        UCVMATH::vec_t ucvmeanv = UCVMATH::vec_new(4);

        for (int i = 0; i < 4; i++)
        {
            ucvmeanv.v[i] = meanArray[i];
        }

        vtkm::IdComponent numSamples = 1000;
        vtkm::Id numCrossings = 0;
        // this can be adapted to 3d case

        UCVMATH::mat ucvcov4by4 = UCVMATH::matrix_new(4, 4);
        int covindex = 0;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                // use the elements at the top half
                // printf("%f ", cov_matrix[covindex]);
                ucvcov4by4->v[p][q] = cov_matrix[covindex];
                if (p != q)
                {
                    // assign value to another helf
                    ucvcov4by4->v[q][p] = ucvcov4by4->v[p][q];
                }
                covindex++;
            }
        }
        //printf("\nshow cov matrix for ucv_matrix func\n");
        //UCVMATH::matrix_show(cov4by4);

        // we have set the eigen value as 0 when it is
        // a really small negative value
        UCVMATH::mat A = UCVMATH::eigen_vector_decomposition(ucvcov4by4);
        //printf("\nucv computing transform matrix:\n");
        //UCVMATH::matrix_show(A);
        //std::cout << std::endl;

        UCVMATH::mat sampleU = UCVMATH::eigen_vector_decomposition(4, numSamples);
        //printf("\nucv computing sampling matrix:\n");
        //UCVMATH::matrix_show(sampleU);

        //printf("\nucv ucvmeanv:\n");
        //UCVMATH::vec_show(&ucvmeanv);

        // computing the AU+M for each column
        UCVMATH::mat AUM = UCVMATH::matrix_mul_add_vector(A, sampleU, &ucvmeanv);
        //printf("\nucv computing AUM:\n");
        //UCVMATH::matrix_show(AUM);

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            if ((m_isovalue <= AUM->v[0][n]) && (m_isovalue <= AUM->v[1][n]) && (m_isovalue <= AUM->v[2][n]) && (m_isovalue <= AUM->v[3][n]))
            {
                numCrossings = numCrossings + 0;
            }
            else if ((m_isovalue >= AUM->v[0][n]) && (m_isovalue >= AUM->v[1][n]) && (m_isovalue >= AUM->v[2][n]) && (m_isovalue >= AUM->v[3][n]))
            {
                numCrossings = numCrossings + 0;
            }
            else
            {
                numCrossings = numCrossings + 1;
            }
        }

        // cross probability
        //std::cout << "ucv numCrossings " << numCrossings << std::endl;
        outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);

        UCVMATH::matrix_delete(AUM);
        UCVMATH::matrix_delete(sampleU);
        UCVMATH::matrix_delete(A);
        UCVMATH::matrix_delete(ucvcov4by4);
        UCVMATH::vec_delete(&ucvmeanv);

       /*
        for debuging, comparing with eigen reuslts
       
        Eigen::Vector4d meanVector(4);
        meanVector << meanArray[0], meanArray[1], meanArray[2], meanArray[3];

        Eigen::Matrix4d cov4by4(4, 4);
        cov4by4 << cov_matrix[0], cov_matrix[1], cov_matrix[2], cov_matrix[3],
             cov_matrix[1], cov_matrix[4], cov_matrix[5], cov_matrix[6],
             cov_matrix[2], cov_matrix[5], cov_matrix[7], cov_matrix[8],
             cov_matrix[3], cov_matrix[6], cov_matrix[8], cov_matrix[9];

        numSamples = 100;
        //compute the eigen version for checking
        Eigen::EigenMultivariateNormal<vtkm::Float64> normX_solver(meanVector, cov4by4);
        auto R = normX_solver.samples(numSamples).transpose();
        int eigenNCrossing=0;
        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            // std::cout << R.coeff(n, 0) << " " << R.coeff(n, 1) << " " << R.coeff(n, 2) << " " << R.coeff(n, 3) << std::endl;
            // TODO, how to make it dynamic for 2d and 3d case
            if ((m_isovalue <= R.coeff(n, 0)) && (m_isovalue <= R.coeff(n, 1)) && (m_isovalue <= R.coeff(n, 2)) && (m_isovalue <= R.coeff(n, 3)))
            {
                eigenNCrossing = eigenNCrossing + 0;
            }
            else if ((m_isovalue >= R.coeff(n, 0)) && (m_isovalue >= R.coeff(n, 1)) && (m_isovalue >= R.coeff(n, 2)) && (m_isovalue >= R.coeff(n, 3)))
            {
                eigenNCrossing = eigenNCrossing + 0;
            }
            else
            {
                eigenNCrossing = eigenNCrossing + 1;
            }
        }

        std::cout << "eigenNCrossing " << eigenNCrossing << std::endl;
        */
    }

    // how to get this vtkm::Vec<double, 15> in an more efficient way
    VTKM_EXEC vtkm::Float64 find_mean(const vtkm::Vec<double, 15> &arr) const
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

    VTKM_EXEC double find_covariance(const vtkm::Vec<double, 15> &arr1, const vtkm::Vec<double, 15> &arr2,
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
        return sum / (1.0 * (arraySize - 1));
    }

private:
    double m_isovalue;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h