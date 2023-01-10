#ifndef UCV_MULTIVARIANT_GAUSSIAN_h
#define UCV_ENTROPY_MULTIVARIANT_GAUSSIAN_h

#include <vector>
#include <vtkm/worklet/WorkletReduceByKey.h>
#include <Eigen/Dense>

class MultivariantGaussian : public vtkm::worklet::WorkletReduceByKey
{
public:
    MultivariantGaussian(double isovalue, vtkm::Id hixelBlockDim)
        : m_isovalue(isovalue),
          m_hixelBlockDim(hixelBlockDim){};

    using ControlSignature = void(KeysIn, ValuesIn, ReducedValuesOut);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;
    template <typename KeyType, typename OriginalValuesType, typename OutputType>
    VTKM_EXEC void operator()(KeyType key,
                              const OriginalValuesType &originalValues, OutputType &outputValue) const

    {

        vtkm::IdComponent numComponents = originalValues.GetNumberOfComponents();

        if (numComponents > m_hixelBlockDim * m_hixelBlockDim)
        {
            throw std::runtime_error("only support 2d case for current MultivariantGaussian implementation");
        }

        // std::cout << "numComponents " << numComponents << std::endl;
        if (numComponents != m_hixelBlockDim * m_hixelBlockDim)
        {
            // this is not whole hixel blocks
            outputValue = 0;
            return;
        }

        vtkm::IdComponent xid, yid, xidnew, yidnew, vertexid;
        // for vertex in each cell
        // each vertex contains hixelBlockDim*hixelBlockDim original data
        std::vector<std::vector<vtkm::FloatDefault>> rawData(4, std::vector<vtkm::FloatDefault>());

        // here we use the small block
        // assuming the actual hixel is 4*4, in order to compute variance
        // we actually have 8*8 data here
        // m_hixelBlockDim is supposed to be divided by 2;
        vtkm::Id actualHixel = m_hixelBlockDim / 2;

        for (vtkm::IdComponent index = 0; index < numComponents; index++)
        {
            // be careful about when to use large hixel block size
            // and when to use small hixel data block size (actual data block size)
            xid = index % m_hixelBlockDim;
            yid = index / m_hixelBlockDim;

            xidnew = xid / actualHixel;
            yidnew = yid / actualHixel;

            // the index for 2*2 grid
            vertexid = yidnew * 2 + xidnew;

            //std::cout << "index " << index << " xid " << xid << " yid " << yid
            //          << " xidnew " << xidnew << " yidnew " << yidnew << " vertexid " << vertexid << std::endl;

            rawData[vertexid].push_back(originalValues[index]);
        }

        // The operations here are same with
        // https://github.com/MengjiaoH/Probabilistic-Marching-Cubes-C-/blob/main/covariance.h
        // compute u00 u10 u01 u11
        std::vector<vtkm::FloatDefault> meanArray(4, 0);
        meanArray[0] = find_mean(rawData[0]);
        meanArray[1] = find_mean(rawData[1]);
        meanArray[2] = find_mean(rawData[2]);
        meanArray[3] = find_mean(rawData[3]);

        // compute cov 4*4
        std::vector<vtkm::FloatDefault> cov_matrix;
        for (int p = 0; p < 4; ++p)
        {
            for (int q = p; q < 4; ++q)
            {
                if (p == q)
                {
                    cov_matrix.push_back(1.0);
                    continue;
                }

                float cov = find_covariance(rawData[p], rawData[q], meanArray[p], meanArray[q]);
                cov_matrix.push_back(cov);
            }
        }

        // generate sample

        // generate mean and cov matrix
        Eigen::Vector4d cov4by4;
        cov4by4 << cov_matrix[0], cov_matrix[1], cov_matrix[2], cov_matrix[3],
            cov_matrix[1], cov_matrix[4], cov_matrix[5], cov_matrix[6],
            cov_matrix[2], cov_matrix[5], cov_matrix[7], cov_matrix[8],
            cov_matrix[3], cov_matrix[6], cov_matrix[8], cov_matrix[9];
        // sample the results from the distribution function and compute the cross probability

        // how to check the accuracy of cov4by4?
    }

    vtkm::FloatDefault find_mean(std::vector<vtkm::FloatDefault> &arr) const
    {
        vtkm::FloatDefault sum = std::accumulate(arr.begin(), arr.end(), 0.0f);
        vtkm::FloatDefault mean = (vtkm::FloatDefault)sum / (vtkm::FloatDefault)(arr.size());
        return mean;
    }

    vtkm::FloatDefault find_covariance(std::vector<vtkm::FloatDefault> &arr1, std::vector<vtkm::FloatDefault> &arr2,
                                       vtkm::FloatDefault &mean1, vtkm::FloatDefault &mean2) const
    {
        if (arr1.size() != arr2.size())
        {
            throw std::runtime_error("failed to compute find_covariance, the array size should be equal with each other");
        }
        vtkm::Id arraySize = arr1.size();
        vtkm::FloatDefault sum = 0;
        for (int i = 0; i < arraySize; i++)
            sum = sum + (arr1[i] - mean1) * (arr2[i] - mean2);
        return sum / (vtkm::FloatDefault)(arraySize - 1);
    }

private:
    double m_isovalue;
    vtkm::Id m_hixelBlockDim;
};

#endif // UCV_ENTROPY_MULTIVARIANT_GAUSSIAN_h