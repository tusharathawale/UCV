#ifndef UCV_CRITICAL_POINT_h
#define UCV_CRITICAL_POINT_h

#include <vtkm/worklet/WorkletPointNeighborhood.h>

struct CriticalPointWorklet : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    CriticalPointWorklet(){};

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
        printf("workIndex is %d\n", WorkIndex);
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

        printf("check input a1 %f b1 %f a2 %f b2 %f a3 %f b3 %f a4 %f b4 %f a5 %f b5 %f\n", a1, b1, a2, b2, a3, b3, a4, b4, a5, b5);

        // compute bmin
        vtkm::FloatDefault bMin = vtkm::Min(b1, vtkm::Min(b2, vtkm::Min(b3, vtkm::Min(b4, b5))));
        printf("bmin %f\n", bMin);

        if (bMin <= a1)
        {
            minProb = 0;
            return;
        }

        // startPointList = [a1, a2, a3, a4, a5]
        // order = np.argsort(startPointList)
        // interval, first is value second is actual index from 0 to 4
        vtkm::Vec<vtkm::Pair<vtkm::FloatDefault, vtkm::FloatDefault>, 5> interval;
        interval[0] = {a1, b1};
        interval[1] = {a2, b2};
        interval[2] = {a3, b3};
        interval[3] = {a4, b4};
        interval[4] = {a5, b5};
        // vtkm::Vec<vtkm::Id, 5> sotedIndex = ArgSort<5>(interval);
        ArgSort<5>(interval);
        printf("sorted a [%f %f %f %f %f]\n", interval[0].first, interval[1].first, interval[2].first, interval[3].first, interval[4].first);

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
        printf("---endInterval %d\n", tartgetIndex);






    }
    // ascending
    template <vtkm::Id Size>
    VTKM_EXEC inline void ArgSort(vtkm::Vec<vtkm::Pair< vtkm::FloatDefault, vtkm::FloatDefault>, Size> &interval) const
    {
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size - i - 1; j++)
            {
                // compare element i and j
                if (interval[j].first > interval[j + 1].first)
                {
                    // change
                    auto temp = interval[j];
                    interval[j] = interval[j + 1];
                    interval[j + 1] = temp;
                }
            }
        }

        // vtkm::Vec<vtkm::Id, Size> result;

        // for (int i = 0; i < Size; i++)
        // {
        //     result[i] = interval[i].first;
        // }
        return ;
    }

private:
    int m_neighborhoodSize = 1;
};

#endif // UCV_CRITICAL_POINT_h