#include <vtkm/cont/Initialize.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/worklet/WorkletPointNeighborhood.h>
#include <vtkm/cont/ArrayHandle.h>

struct ComputeMoments3D : public vtkm::worklet::WorkletPointNeighborhood
{
public:
    ComputeMoments3D(const vtkm::Vec3f &_spacing, vtkm::Float64 _radius, int _p, int _q, int _r)
        : RadiusDiscrete(vtkm::IdComponent(_radius / (_spacing[0] - 1e-10)),
                         vtkm::IdComponent(_radius / (_spacing[1] - 1e-10)),
                         vtkm::IdComponent(_radius / (_spacing[2] - 1e-10))),
          SpacingProduct(vtkm::ReduceProduct(_spacing)), p(_p), q(_q), r(_r)
    {
        assert(_spacing[0] > 1e-10);
        assert(_spacing[1] > 1e-10);
        assert(_spacing[2] > 1e-10);

        assert(_p >= 0);
        assert(_q >= 0);
        assert(_r >= 0);
    }

    using ControlSignature = void(CellSetIn, FieldInNeighborhood, FieldOut);

    using ExecutionSignature = void(_2, Boundary, _3);

    template <typename NeighIn, typename T>
    VTKM_EXEC void operator()(const NeighIn &image,
                              const vtkm::exec::BoundaryState &boundary,
                              T &moment) const
    {
        // TODO: type safety and numerical precision
        auto sum = vtkm::TypeTraits<T>::ZeroInitialization();

        // Clamp the radius to the dataset bounds (discard out-of-bounds points).
        const auto minRadius = boundary.ClampNeighborIndex(-this->RadiusDiscrete);
        const auto maxRadius = boundary.ClampNeighborIndex(this->RadiusDiscrete);
        std::cout << "k dim " << minRadius[2] << "," << maxRadius[2] << std::endl;
        std::cout << "j dim " << minRadius[1] << "," << maxRadius[1] << std::endl;
        std::cout << "i dim " << minRadius[0] << "," << maxRadius[0] << std::endl;

        vtkm::Vec3f_64 radius;
        for (vtkm::IdComponent k = minRadius[2]; k <= maxRadius[2]; ++k)
        {
            if (k > -this->RadiusDiscrete[2] && boundary.IJK[2] + k == 0)
            { // Don't double count samples that exist on other nodes:
                continue;
            }
            radius[2] = k * 1. / this->RadiusDiscrete[2];

            for (vtkm::IdComponent j = minRadius[1]; j <= maxRadius[1]; ++j)
            {
                if (j > -this->RadiusDiscrete[1] && boundary.IJK[1] + j == 0)
                { // Don't double count samples that exist on other nodes:
                    continue;
                }
                radius[1] = j * 1. / this->RadiusDiscrete[1];

                for (vtkm::IdComponent i = minRadius[0]; i <= maxRadius[0]; ++i)
                {
                    if (i > -this->RadiusDiscrete[0] && boundary.IJK[0] + i == 0)
                    { // Don't double count samples that exist on other nodes:
                        continue;
                    }
                    radius[0] = i * 1. / this->RadiusDiscrete[0];

                    if (vtkm::Dot(radius, radius) <= 1)
                    {
                        sum += static_cast<T>(vtkm::Pow(radius[0], p) * vtkm::Pow(radius[1], q) *
                                              vtkm::Pow(radius[2], r) * image.Get(i, j, k));
                    }
                }
            }
        }

        moment = T(sum * this->SpacingProduct);
    }

private:
    vtkm::Vec3i_32 RadiusDiscrete;
    const vtkm::Float64 SpacingProduct;
    const int p;
    const int q;
    const int r;
};

using SupportedTypes = vtkm::List<vtkm::Float32,
                                  vtkm::Float64,
                                  vtkm::Int8,
                                  vtkm::UInt8,
                                  vtkm::Int16,
                                  vtkm::UInt16,
                                  vtkm::Int32,
                                  vtkm::UInt32>;

int main(int argc, char *argv[])
{
    // init the vtkm (set the backend and log level here)
    vtkm::cont::Initialize(argc, argv);

    if (argc != 2)
    {
        std::cout << "executable <filename>" << std::endl;
        exit(0);
    }

    std::string fileName = argv[1];
    // load the dataset (beetles data set, structured one)
    // TODO, the data set can be distributed between different ranks

    // create the vtkm data set from the loaded data
    vtkm::io::VTKDataSetReader reader(fileName);
    vtkm::cont::DataSet inData = reader.ReadDataSet();

    // check the property of the data
    // the raw_data for testing is 832*832*494
    inData.PrintSummary(std::cout);

    // TODO try to use the point neighborhood
    // skip particular on if it is covered by previous points
    // refer to this:
    // https://gitlab.kitware.com/vtk/vtk-m/-/blob/release-1.9/vtkm/filter/image_processing/worklet/ComputeMoments.h

    /* try the neighborhood worklet */

    std::cout << "get cell set" << std::endl;
    const vtkm::cont::UnknownCellSet &inputCellSet = inData.GetCellSet();
    std::cout << "get field" << std::endl;

    auto field = inData.GetField("ground_truth");

    vtkm::cont::UnknownArrayHandle outArray;

    using WorkletType = ComputeMoments3D;
    using DispatcherType = vtkm::worklet::DispatcherPointNeighborhood<WorkletType>;

    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> result;
    std::cout << "dispatcher and invoke" << std::endl;

    auto resolveType = [&](const auto &concrete)
    {
        DispatcherType dispatcher(WorkletType{{0.1, 0.1, 0.1}, 2, 1, 1, 1});
        dispatcher.Invoke(inputCellSet, concrete, result);
    };

    field.GetData().CastAndCallForTypesWithFloatFallback<SupportedTypes, VTKM_DEFAULT_STORAGE_LIST>(
        resolveType);

    return 0;
}
