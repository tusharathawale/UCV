// this is a temporary file for another purpose

#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleCompositeVector.h>
#include <vtkm/cont/DataSet.h>

#include <vector>

int main(int argc, char *argv[])
{
    // create a data set with multiple fields
    vtkm::cont::DataSet dataSet;
    std::vector<vtkm::FloatDefault> pointarray;
    for (vtkm::Id i = 0; i < 100; ++i)
    {
        pointarray.push_back(i * 0.1);
    }

    dataSet.AddPointField("field1", pointarray);
    dataSet.AddPointField("field2", pointarray);
    dataSet.AddPointField("field3", pointarray);

    // composite the field into the composite vector

    vtkm::cont::UnknownArrayHandle results;

    auto resolveType3d = [&](const auto &field1)
    {
        // assume three field have the same type
        using T = typename std::decay_t<decltype(field1)>::ValueType;
        vtkm::cont::ArrayHandle<T> field2;
        vtkm::cont::ArrayCopyShallowIfPossible(dataSet.GetField("field2").GetData(), field2);
        vtkm::cont::ArrayHandle<T> field3;
        vtkm::cont::ArrayCopyShallowIfPossible(dataSet.GetField("field3").GetData(), field3);

        auto composite = vtkm::cont::make_ArrayHandleCompositeVector(field1, field2, field3);

        results = composite;
    };

    const auto &inField1 = dataSet.GetField("field1");

    inField1.GetData().CastAndCallForTypesWithFloatFallback<vtkm::TypeListFieldScalar, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType3d);

    vtkm::cont::Field outField3d("composite3d",
                                 inField1.GetAssociation(),
                                 results);
    dataSet.AddField(outField3d);

    auto resolveType2d = [&](const auto &field1)
    {
        // assume three field have the same type
        using T = typename std::decay_t<decltype(field1)>::ValueType;
        vtkm::cont::ArrayHandle<T> field2;
        vtkm::cont::ArrayCopyShallowIfPossible(dataSet.GetField("field2").GetData(), field2);

        auto composite = vtkm::cont::make_ArrayHandleCompositeVector(field1, field2);

        results = composite;
    };

    inField1.GetData().CastAndCallForTypesWithFloatFallback<vtkm::TypeListFieldScalar, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType2d);

    vtkm::cont::Field outField("composite2d",
                               inField1.GetAssociation(),
                               results);
    dataSet.AddField(outField);

    // check the output data
    dataSet.PrintSummary(std::cout);

    return 0;
}