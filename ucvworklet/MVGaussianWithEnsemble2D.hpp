#ifndef UCV_MULTIVARIANT_GAUSSIAN2D_h
#define UCV_MULTIVARIANT_GAUSSIAN2D_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>
#include <Eigen/Dense>
#include "./eigenmvn.h"
// this worklet is for the input data that put the different data in a separate array
// for the wind data here https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/tree/main/datasets/txt_files/wind_pressure_200
// there are 15 numbers (ensemble extraction) each data is put in a different file

class MVGaussianWithEnsemble2D : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble2D(double isovalue)
        : m_isovalue(isovalue){};

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
        OutCellFieldType &outCellFieldCProb,
        vtkm::Id workIndex) const
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

        std::vector<vtkm::Float64> meanArray(4, 0);

        // derive type
        meanArray[0] = find_mean(inPointFieldVecEnsemble[updateIndex4(0)]);
        meanArray[1] = find_mean(inPointFieldVecEnsemble[updateIndex4(1)]);
        meanArray[2] = find_mean(inPointFieldVecEnsemble[updateIndex4(2)]);
        meanArray[3] = find_mean(inPointFieldVecEnsemble[updateIndex4(3)]);

        //if (workIndex == 0)
        //{
        //    std::cout << "mean " << meanArray[0] << " " << meanArray[1] << " " << meanArray[2] << " " << meanArray[3] << std::endl;
        //}

        std::vector<float> cov_matrix;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                int updatep = updateIndex4(p);
                int updateq = updateIndex4(q);
                float cov = find_covariance(inPointFieldVecEnsemble[updatep], inPointFieldVecEnsemble[updateq], meanArray[p], meanArray[q]);
                cov_matrix.push_back(cov);
                if (workIndex == 0)
                {
                    std::cout << cov << " " << std::endl;
                }
            }
        }

        // generate sample

        Eigen::Vector4d meanVector(4);
        meanVector << meanArray[0], meanArray[1], meanArray[2], meanArray[3];

        // for (int i = 0; i < 4; i++)
        //{
        //     meanVector(i) = meanArray[i];
        // }

        // std::cout << "cov_matrix size " << cov_matrix.size() << std::endl;

        // generate mean and cov matrix
        Eigen::Matrix4d cov4by4(4, 4);
        cov4by4 << cov_matrix[0], cov_matrix[1], cov_matrix[2], cov_matrix[3],
            cov_matrix[1], cov_matrix[4], cov_matrix[5], cov_matrix[6],
            cov_matrix[2], cov_matrix[5], cov_matrix[7], cov_matrix[8],
            cov_matrix[3], cov_matrix[6], cov_matrix[8], cov_matrix[9];
        // this can be adapted to 3d case
        /*

        int covindex = 0;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                // use the elements at the top half
                cov4by4(p, q) = cov_matrix[covindex];
                if (p != q)
                {
                    // assign value to another helf
                    cov4by4(q, p) = cov4by4(p, q);
                }
                covindex++;
            }
        }
        */

        // sample the results from the distribution function and compute the cross probability
        vtkm::IdComponent numSamples = 1000;
        Eigen::EigenMultivariateNormal<vtkm::Float64> normX_solver(meanVector, cov4by4);

        auto R = normX_solver.samples(numSamples).transpose();

        vtkm::Id numCrossings = 0;

        for (vtkm::Id n = 0; n < numSamples; ++n)
        {
            //if (workIndex == 0 && n==0)
            //{
            //    std::cout << "check r" << std::endl;
            //    std::cout << R.coeff(0, 0) << " " << R.coeff(0, 1) << " " << R.coeff(0, 2) << " " << R.coeff(0, 3) << std::endl;
            //}
            //
            // TODO, how to make it dynamic for 2d and 3d case
            if ((m_isovalue <= R.coeff(n, 0)) && (m_isovalue <= R.coeff(n, 1)) && (m_isovalue <= R.coeff(n, 2)) && (m_isovalue <= R.coeff(n, 3)))
            {
                numCrossings = numCrossings + 0;
            }
            else if ((m_isovalue >= R.coeff(n, 0)) && (m_isovalue >= R.coeff(n, 1)) && (m_isovalue >= R.coeff(n, 2)) && (m_isovalue >= R.coeff(n, 3)))
            {
                numCrossings = numCrossings + 0;
            }
            else
            {
                numCrossings = numCrossings + 1;
            }
        }

        // cross probability
        outCellFieldCProb = (double)numCrossings / numSamples;
    }

    int updateIndex4(int index) const
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

    // how to get this vtkm::Vec<double, 15> in an more efficient way
    vtkm::Float64 find_mean(const vtkm::Vec<double, 15> &arr) const
    {
        vtkm::Float64 sum = 0;
        int num = 15;

        for (vtkm::Id i = 0; i < 15; i++)
        {
            sum = sum + arr[i];
        }
        vtkm::Float64 mean = (double)sum / (double)num;
        return mean;
    }

    double find_covariance(const vtkm::Vec<double, 15> &arr1, const vtkm::Vec<double, 15> &arr2,
                           double &mean1, double &mean2) const
    {
        if (arr1.GetNumberOfComponents() != arr2.GetNumberOfComponents())
        {
            // cuda does not support exception
            printf("error, failed to compute find_covariance, the array size should be equal with each other\n");
            return 0;
        }
        vtkm::Id arraySize = arr1.GetNumberOfComponents();
        double sum = 0.0;
        for (int i = 0; i < arraySize; i++)
            sum = sum + (arr1[i] - mean1) * (arr2[i] - mean2);
        return sum / (double)(arraySize - 1);
    }

private:
    double m_isovalue;
};

#endif // UCV_MULTIVARIANT_GAUSSIAN2D_h