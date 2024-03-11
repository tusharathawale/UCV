#ifndef UCV_CRITICAL_POINT_AVOID_OVERFLOW_h
#define UCV_CRITICAL_POINT_AVOID_OVERFLOW_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>

#ifdef USE_LOG
#define LOG(x) x
#else
#define LOG(x)
#endif

struct CriticalPointWorkletAvoidOverflow : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    CriticalPointWorkletAvoidOverflow(){};

    using ControlSignature = void(CellSetIn, FieldInNeighborhood, FieldInNeighborhood, FieldOut);

    using ExecutionSignature = void(_2, _3, _4, Boundary, WorkIndex);

    template <typename InPointField, typename OutPointField>
    VTKM_EXEC void operator()(const InPointField &minValue,
                              const InPointField &maxValue,
                              OutPointField &minProb,
                              const vtkm::exec::BoundaryState &boundary,
                              vtkm::Id WorkIndex) const
    {
        // resluts is the coordinates of three dims
        auto minIndices = boundary.MinNeighborIndices(this->m_neighborhoodSize);
        auto maxIndices = boundary.MaxNeighborIndices(this->m_neighborhoodSize);

        // minIndices is supposed to be -1
        // maxIndices is supposed to be 1
        // if (WorkIndex == 0)
        //{
        // debug
        LOG(printf("workIndex is %d\n", WorkIndex));
        // printf("min index %d %d %d\n", minIndices[0], minIndices[1], minIndices[2]);
        // printf("max index %d %d %d\n", maxIndices[0], maxIndices[1], maxIndices[2]);

        // filter out the element in the boundry
        // if the element is at the boundry places its min prob is 0
        if ((maxIndices[0] - minIndices[0] < 2) || (maxIndices[1] - minIndices[1] < 2))
        {
            // if x and y is at the boundry, do not consider it
            minProb = 0;
            return;
        }
        // for testing
        // minProb=1;
        // get a1-a5 b1-b5
        //  a1 b1 self
        //  a2 b2 i+1 j
        //  a3 b3 i-1 j
        //  a4 b4 i,j+1
        //  a5 b5 i,j-1
        // in vtkm the i and j is inverted compared with python

        vtkm::FloatDefault a1 = minValue.Get(0, 0, 0);
        vtkm::FloatDefault b1 = maxValue.Get(0, 0, 0);

        vtkm::FloatDefault a2 = minValue.Get(0, 1, 0);
        vtkm::FloatDefault b2 = maxValue.Get(0, 1, 0);

        vtkm::FloatDefault a3 = minValue.Get(0, -1, 0);
        vtkm::FloatDefault b3 = maxValue.Get(0, -1, 0);

        vtkm::FloatDefault a4 = minValue.Get(1, 0, 0);
        vtkm::FloatDefault b4 = maxValue.Get(1, 0, 0);

        vtkm::FloatDefault a5 = minValue.Get(-1, 0, 0);
        vtkm::FloatDefault b5 = maxValue.Get(-1, 0, 0);

        LOG(printf("check input a1 %f b1 %f a2 %f b2 %f a3 %f b3 %f a4 %f b4 %f a5 %f b5 %f\n", a1, b1, a2, b2, a3, b3, a4, b4, a5, b5));
       
        //do the preprocessing to avoid the data overflow
        vtkm::FloatDefault a1N = (a1-a1)*this->m_ScaleNum;
        vtkm::FloatDefault b1N = (b1-a1)*this->m_ScaleNum;

        vtkm::FloatDefault a2N = (a2-a1)*this->m_ScaleNum;
        vtkm::FloatDefault b2N = (b2-a1)*this->m_ScaleNum;

        vtkm::FloatDefault a3N = (a3-a1)*this->m_ScaleNum;
        vtkm::FloatDefault b3N = (b3-a1)*this->m_ScaleNum;

        vtkm::FloatDefault a4N = (a4-a1)*this->m_ScaleNum;
        vtkm::FloatDefault b4N = (b4-a1)*this->m_ScaleNum;
        
        vtkm::FloatDefault a5N = (a5-a1)*this->m_ScaleNum;
        vtkm::FloatDefault b5N = (b5-a1)*this->m_ScaleNum;

        LOG(printf("check transformed input a1 %f b1 %f a2 %f b2 %f a3 %f b3 %f a4 %f b4 %f a5 %f b5 %f\n", a1N, b1N, a2N, b2N, a3N, b3N, a4N, b4N, a5N, b5N));

        // compute bmin
        vtkm::FloatDefault bMin = vtkm::Min(b1N, vtkm::Min(b2N, vtkm::Min(b3N, vtkm::Min(b4N, b5N))));
        LOG(printf("bmin %f\n", bMin));

        if (bMin <= a1N)
        {
            minProb = 0;
            return;
        }

        // startPointList = [a1, a2, a3, a4, a5]
        // order = np.argsort(startPointList)
        // interval, first is value second is actual index from 0 to 4
        vtkm::Vec<vtkm::Pair<vtkm::Float64, vtkm::Float64>, 5> interval;
        interval[0] = {a1N, b1N};
        interval[1] = {a2N, b2N};
        interval[2] = {a3N, b3N};
        interval[3] = {a4N, b4N};
        interval[4] = {a5N, b5N};
        // vtkm::Vec<vtkm::Id, 5> sotedIndex = ArgSort<5>(interval);
        ArgSort<5>(interval);
        LOG(printf("sorted a [%f %f %f %f %f]\n", interval[0].first, interval[1].first, interval[2].first, interval[3].first, interval[4].first));

        // find interval contain bMin
        // the interval vector is sorted now
        vtkm::Id tartgetIndex;
        for (tartgetIndex = 4; tartgetIndex > 0; tartgetIndex--)
        {
            auto intervalFromEnd = interval[tartgetIndex];
            if ((bMin >= intervalFromEnd.first) && (bMin <= intervalFromEnd.second))
            {
                break;
            }
        }
        LOG(printf("---endInterval %d\n", tartgetIndex));

        // create x1 limit, init it as negative value
        vtkm::Vec<vtkm::Float64, 6> x1Limit(vtkm::Nan64());
        vtkm::Id index;
        for (index = 0; index < tartgetIndex + 1; index++)
        {
            x1Limit[index] = interval[index].first;
        }
        // add bMin as the last element
        x1Limit[index] = bMin;

        // find the a1
        vtkm::Id indexa1 = -1;
        for (vtkm::Id i = 0; i < 6; i++)
        {
            if (vtkm::Abs(x1Limit[i] - a1N) < 0.0000001)
            {
                indexa1 = i;
                break;
            }
        }

        if (indexa1 == -1)
        {
            // run time error
            printf("error, the indexa1 is not supposed to be -1");
            minProb = 0;
            return;
        }

        vtkm::FloatDefault w1 = b1N - a1N;
        // call superOptimizedCase
        minProb = SuperOptimizedCase(indexa1, x1Limit, interval, w1);
        return;
    }
    // ascending
    template <vtkm::Id Size>
    VTKM_EXEC inline void ArgSort(vtkm::Vec<vtkm::Pair<vtkm::Float64, vtkm::Float64>, Size> &interval) const
    {
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size - i - 1; j++)
            {
                // compare element i and j
                if (interval[j].first > interval[j + 1].first)
                {
                    // swap
                    auto temp = interval[j];
                    interval[j] = interval[j + 1];
                    interval[j + 1] = temp;
                }
            }
        }
        return;
    }
    //        minimaProb = superOptimizedCase(indexOfa1,x1Limits, sortedI, w1)
    VTKM_EXEC inline vtkm::Float64 SuperOptimizedCase(vtkm::Id indexOfa1,
                                                      vtkm::Vec<vtkm::Float64, 6> x1Limits,
                                                      vtkm::Vec<vtkm::Pair<vtkm::Float64, vtkm::Float64>, 5> interval,
                                                      vtkm::Float64 w1) const
    {
        LOG(printf("---debug SuperOptimizedCase indexOfa1 %d w1 %f \n", indexOfa1, w1));
        LOG(printf("---debug x1Limits\n"));
        for (int i = 0; i < 6; i++)
        {
            LOG(printf("%f ", x1Limits[i]));
        }
        LOG(printf("\n"));
        LOG(printf("---debug interval\n"));
        for (int i = 0; i < 5; i++)
        {
            LOG(printf("[%f %f]", interval[i].first, interval[i].second));
        }
        LOG(printf("\n"));
        vtkm::Float64 minProb = 0;
        // upLimits = [0,None,None,None,None,None]
        // normalizerWidths = [0, w1, 1,1,1,1]
        vtkm::Vec<vtkm::Float64, 6> upLimits({0, vtkm::Nan64(), vtkm::Nan64(), vtkm::Nan64(), vtkm::Nan64(), vtkm::Nan64()});
        vtkm::Vec<vtkm::Float64, 6> normalizerWidths({0, w1, 1, 1, 1, 1});

        for (vtkm::Id k = 0; k < indexOfa1; k++)
        {
            auto tempInterval = interval[k];
            upLimits[k + 2] = tempInterval.second;
            normalizerWidths[k + 2] = tempInterval.second - tempInterval.first;
        }

        // computing cases before indexOfa1
        minProb = minProb + ComputeIntegral(x1Limits[indexOfa1], x1Limits[indexOfa1 + 1], upLimits[2], upLimits[3], upLimits[4], upLimits[5], normalizerWidths[1], normalizerWidths[2], normalizerWidths[3], normalizerWidths[4], normalizerWidths[5]);
        LOG(printf("----debug minProb first case %f\n", minProb));
        // computing cases after indexOfa1
        for (vtkm::Id i = indexOfa1 + 1; i < 5; i++)
        {
            // break if x1Limits goes to nan
            if (x1Limits[i] == vtkm::Nan64())
            {
                break;
            }
            auto tempInterval = interval[i];
            upLimits[i + 1] = tempInterval.second;
            normalizerWidths[i + 1] = tempInterval.second - tempInterval.first;
            minProb = minProb + ComputeIntegral(x1Limits[i], x1Limits[i + 1], upLimits[2], upLimits[3], upLimits[4], upLimits[5], normalizerWidths[1], normalizerWidths[2], normalizerWidths[3], normalizerWidths[4], normalizerWidths[5]);
        }
        return minProb;
    }

    VTKM_EXEC inline vtkm::Float64 ComputeIntegral(vtkm::Float64 l,
                                                   vtkm::Float64 h,
                                                   vtkm::Float64 h2,
                                                   vtkm::Float64 h3,
                                                   vtkm::Float64 h4,
                                                   vtkm::Float64 h5,
                                                   vtkm::Float64 n1,
                                                   vtkm::Float64 n2,
                                                   vtkm::Float64 n3,
                                                   vtkm::Float64 n4,
                                                   vtkm::Float64 n5) const
    {
        LOG(printf("---debug ComputeIntegral input %f %f %f %f %f %f %f %f %f %f %f\n", l, h, h2, h3, h4, h5, n1, n2, n3, n4, n5));
        vtkm::Float64 intUp = 0;
        vtkm::Float64 intDown = 0;

        bool ln = isnan(l);
        bool hn = isnan(h);
        bool h2n = isnan(h2);
        bool h3n = isnan(h3);
        bool h4n = isnan(h4);
        bool h5n = isnan(h5);
        bool n1n = isnan(n1);
        bool n2n = isnan(n2);
        bool n3n = isnan(n3);
        bool n4n = isnan(n4);
        bool n5n = isnan(n5);

        // printf("check if nan %d %d %d\n",isnan(h3), isnan(h4),isnan(h5));
        vtkm::Float64 normalizingFactor = 1.0 / (n1 * n2 * n3 * n4 * n5);

        if ((!ln) && (!hn) && (h2n ) && (h3n) && (h4n) && (h5n))
        {
            LOG(printf("---c1\n"));
            intUp = normalizingFactor * h;
            intDown = normalizingFactor * l;
        }

        if ((!ln) && (!hn ) && (!h2n) && (h3n) && (h4n) && (h5n))
        {
            LOG(printf("---c2\n"));
            intUp = normalizingFactor * (h2 * h - h * h / 2);
            intDown = normalizingFactor * (h2 * l - l * l / 2);
        }

        if ((!ln) && (!hn ) && (!h2n ) && (!h3n) && (h4n) && (h5n))
        {
            LOG(printf("---c3\n"));
            intUp = normalizingFactor * (h3 * h2 * h - (h3 + h2) * h * h / 2 + h * h * h / 3);
            intDown = normalizingFactor * (h3 * h2 * l - (h3 + h2) * l * l / 2 + l * l * l / 3);
        }

        if ((!ln ) && (!hn ) && (!h2n ) && (!h3n) && (!h4n) && (h5n))
        {
            LOG(printf("---c4\n"));
            intUp = normalizingFactor * (h4 * h3 * h2 * h - (h2 * h3 + h2 * h4 + h3 * h4) * (h * h / 2) + (h2 + h3 + h4) * (h * h * h / 3) - h * h * h * h / 4);
            intDown = normalizingFactor * (h4 * h3 * h2 * l - (h2 * h3 + h2 * h4 + h3 * h4) * (l * l / 2) + (h2 + h3 + h4) * (l * l * l / 3) - l * l * l * l / 4);
        }

        if ((!ln) && (!hn) && (!h2n) && (!h3n ) && (!h4n) && (!h5n))
        {
            LOG(printf("---c5\n"));
            intUp = normalizingFactor * (h5 * h4 * h3 * h2 * h - (h2 * h3 * h4 + h2 * h3 * h5 + h2 * h4 * h5 + h3 * h4 * h5) * (h * h / 2) + (h2 * h3 + h2 * h4 + h2 * h5 + h3 * h4 + h3 * h5 + h4 * h5) * (h * h * h / 3) - (h2 + h3 + h4 + h5) * (h * h * h * h / 4) + h * h * h * h * h / 5);
            intDown = normalizingFactor * (h5 * h4 * h3 * h2 * l - (h2 * h3 * h4 + h2 * h3 * h5 + h2 * h4 * h5 + h3 * h4 * h5) * (l * l / 2) + (h2 * h3 + h2 * h4 + h2 * h5 + h3 * h4 + h3 * h5 + h4 * h5) * (l * l * l / 3) - (h2 + h3 + h4 + h5) * (l * l * l * l / 4) + l * l * l * l * l / 5);
        }

        return (intUp - intDown);
    }

private:
    int m_neighborhoodSize = 1;
    double m_ScaleNum = 1000;
};

#endif // UCV_CRITICAL_POINT_h