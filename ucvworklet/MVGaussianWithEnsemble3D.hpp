#ifndef UCV_MULTIVARIANT_GAUSSIAN3D_h
#define UCV_MULTIVARIANT_GAUSSIAN3D_h

#include <vtkm/worklet/WorkletMapTopology.h>
#include <cmath>
#include <Eigen/Dense>
#include "./eigenmvn.h"

class MVGaussianWithEnsemble3D : public vtkm::worklet::WorkletVisitCellsWithPoints
{
public:
    MVGaussianWithEnsemble3D(double isovalue, int numSamples)
        : m_isovalue(isovalue),m_numSamples(numSamples){};

    using ControlSignature = void(CellSetIn,
                                  FieldInPoint,
                                  FieldInPoint,
                                  FieldOutCell);

    using ExecutionSignature = void(_2, _3, _4);

    // the first parameter is binded with the worklet
    using InputDomain = _1;
    // InPointFieldType should be a vector
    template <typename InPointFieldVecEnsemble, typename InPointFieldVecMean, typename OutCellFieldType>

    VTKM_EXEC void operator()(
        const InPointFieldVecEnsemble &inPointFieldVecEnsemble,
        const InPointFieldVecMean &inMeanArray,
        OutCellFieldType &outCellFieldCProb) const
    {
        // how to process the case where there are multiple variables
        vtkm::IdComponent numVertexies = inPointFieldVecEnsemble.GetNumberOfComponents();

        if (numVertexies != 8)
        {
            //throw std::runtime_error("MVGaussianWithEnsemble3D expects 8 vertecies");
            printf("MVGaussianWithEnsemble3D expects 8 vertecies\n");
            return;
        }

        if (inMeanArray.GetNumberOfComponents() != 8)
        {
            //throw std::runtime_error("inMeanArray in MVGaussianWithEnsemble3D expects 8 vertecies");
            printf("inMeanArray in MVGaussianWithEnsemble3D expects 8 vertecies\n");
            return;
        }

        if (inPointFieldVecEnsemble[0].GetNumberOfComponents() != 64)
        {
            //throw std::runtime_error("only support ensemble size 64 for blockSize equals to 4");
            printf("only support ensemble size 64 for blockSize equals to 4\n");
            return;
        }


        std::vector<double> cov_matrix;
        for (int p = 0; p < numVertexies; ++p)
        {
            for (int q = p; q < numVertexies; ++q)
            {
                float cov = find_covariance(inPointFieldVecEnsemble[p], inPointFieldVecEnsemble[q], inMeanArray[p], inMeanArray[q]);
                cov_matrix.push_back(cov);
            }
        }

        // generate sample

        Eigen::VectorXd meanVector(numVertexies);
        for (int i = 0; i < numVertexies; i++)
        {
            meanVector(i) = inMeanArray[i];
        }

        // generate mean and cov matrix
        Eigen::MatrixXd covMatrix(numVertexies, numVertexies);
        ;
        int covindex = 0;
        for (int p = 0; p < numVertexies; ++p)
        {
            for (int q = p; q < numVertexies; ++q)
            {
                // use the elements at the top half
                covMatrix(p, q) = cov_matrix[covindex];
                if (p != q)
                {
                    // assign value to another helf
                    covMatrix(q, p) = covMatrix(p, q);
                }
                covindex++;
            }
        }
        // sample the results from the distribution function and compute the cross probability
        vtkm::IdComponent numSamples = this->m_numSamples;
        Eigen::EigenMultivariateNormal<vtkm::FloatDefault> normX_solver(meanVector, covMatrix);

        auto R = normX_solver.samples(numSamples).transpose();

        vtkm::Id numCrossings = 0;

        for (int n = 0; n < numSamples; ++n)
        {
            // std::cout << R.coeff(n, 0) << " " << R.coeff(n, 1) << " " << R.coeff(n, 2) << " " << R.coeff(n, 3) << std::endl;

            if ((m_isovalue <= R.coeff(n, 0)) && (m_isovalue <= R.coeff(n, 1)) && (m_isovalue <= R.coeff(n, 2)) && (m_isovalue <= R.coeff(n, 3)) && (m_isovalue <= R.coeff(n, 4)) && (m_isovalue <= R.coeff(n, 5)) && (m_isovalue <= R.coeff(n, 6)) && (m_isovalue <= R.coeff(n, 7)))
            {
                numCrossings = numCrossings + 0;
            }
            else if ((m_isovalue >= R.coeff(n, 0)) && (m_isovalue >= R.coeff(n, 1)) && (m_isovalue >= R.coeff(n, 2)) && (m_isovalue >= R.coeff(n, 3)) && (m_isovalue >= R.coeff(n, 4)) && (m_isovalue >= R.coeff(n, 5)) && (m_isovalue >= R.coeff(n, 6)) && (m_isovalue >= R.coeff(n, 7)))
            {
                numCrossings = numCrossings + 0;
            }
            else
            {
                numCrossings = numCrossings + 1;
            }
        }
        // cross probability
        outCellFieldCProb = (1.0 * numCrossings) / (1.0 * numSamples);
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

    vtkm::FloatDefault find_covariance(const vtkm::Vec<vtkm::FloatDefault, 64> &arr1, const vtkm::Vec<vtkm::FloatDefault, 64> &arr2,
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