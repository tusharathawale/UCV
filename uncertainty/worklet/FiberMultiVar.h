#ifndef vtk_m_worklet_uncertainty_Fiber_h
#define vtk_m_worklet_uncertainty_Fiber_h
#include <random>
#include <iostream>
#include <utility>
#include <vector>
#include <vtkm/worklet/WorkletPointNeighborhood.h>

namespace vtkm
{
namespace worklet
{
namespace detail
{

class VTKM_FILTER_UNCERTAINTY_EXPORT FiberMultiVar : public vtkm::worklet::WorkletPointNeighborhood
{
public:
  // Constructor with initializer list
  FiberMultiVar(const vtkm::Vec<vtkm::Float64, 3>& bottomLeft,
                const vtkm::Vec<vtkm::Float64, 3>& topRight)
    : InputBottomLeft(bottomLeft)
    , InputTopRight(topRight)
  {
  }

  // Input and Output Parameters
  using ControlSignature = void(CellSetIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldIn, FieldOut, FieldOut);
  using ExecutionSignature = void(_2, _3, _4, _5, _6, _7, _8, _9);
  using InputDomain = _1;

  // Template
  template <typename MinX,
            typename MaxX,
            typename MinY,
            typename MaxY,
            typename MinZ,
            typename MaxZ,
            typename OutCellFieldType1,
            typename OutCellFieldType2>
  // Operator
  VTKM_EXEC void operator()(const MinX& EnsembleMinX,
                            const MaxX& EnsembleMaxX,
                            const MinY& EnsembleMinY,
                            const MaxY& EnsembleMaxY,
                            const MinZ& EnsembleMinZ,
                            const MaxZ& EnsembleMaxZ,
                            OutCellFieldType1& MonteCarloProbability,
                            OutCellFieldType2& InteriorProbability) const
  {
    vtkm::FloatDefault minX_user = static_cast<vtkm::FloatDefault>(InputBottomLeft[0]);
    vtkm::FloatDefault minY_user = static_cast<vtkm::FloatDefault>(InputBottomLeft[1]);
    vtkm::FloatDefault minZ_user = static_cast<vtkm::FloatDefault>(InputBottomLeft[2]);
    
    vtkm::FloatDefault maxX_user = static_cast<vtkm::FloatDefault>(InputTopRight[0]);
    vtkm::FloatDefault maxY_user = static_cast<vtkm::FloatDefault>(InputTopRight[1]);
    vtkm::FloatDefault maxZ_user = static_cast<vtkm::FloatDefault>(InputTopRight[2]);

    vtkm::FloatDefault minX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinX);
    vtkm::FloatDefault minY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinY);
    vtkm::FloatDefault minZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMinZ);

    vtkm::FloatDefault maxX_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxX);
    vtkm::FloatDefault maxY_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxY);
    vtkm::FloatDefault maxZ_dataset = static_cast<vtkm::FloatDefault>(EnsembleMaxZ);

    // bottom left
    vtkm::FloatDefault minX_intersection = std::max(minX_user, minX_dataset);
    vtkm::FloatDefault minY_intersection = std::max(minY_user, minY_dataset);
    vtkm::FloatDefault minZ_intersection = std::max(minZ_user, minZ_dataset);

    // top right
    vtkm::FloatDefault maxX_intersection = std::min(maxX_user, maxX_dataset);
    vtkm::FloatDefault maxY_intersection = std::min(maxY_user, maxY_dataset);
    vtkm::FloatDefault maxZ_intersection = std::min(maxZ_user, maxZ_dataset);

    if ((minX_intersection < maxX_intersection) and (minY_intersection < maxY_intersection) and (minZ_intersection < maxZ_intersection))
    {
      vtkm::FloatDefault volume_intersection = (maxX_intersection - minX_intersection) * (maxY_intersection - minY_intersection) * (maxZ_intersection - minZ_intersection);
      vtkm::FloatDefault volume_data = (maxX_dataset - minX_dataset) * (maxY_dataset - minY_dataset) * (maxZ_dataset - minZ_dataset);
      InteriorProbability = volume_intersection / volume_data;
    }
    else
    {
      InteriorProbability = 0.0;
    }

    vtkm::FloatDefault X = 0.0;
    vtkm::FloatDefault Y = 0.0;
    vtkm::FloatDefault Z = 0.0;
    vtkm::IdComponent NumSample = 100;
    vtkm::IdComponent NonZeroCases = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<vtkm::FloatDefault> GenerateX(minX_dataset, maxX_dataset);
    std::uniform_real_distribution<vtkm::FloatDefault> GenerateY(minY_dataset, maxY_dataset);
    std::uniform_real_distribution<vtkm::FloatDefault> GenerateZ(minZ_dataset, maxZ_dataset);
    for (vtkm::IdComponent i = 0; i < NumSample; i++)
    {
      X = GenerateX(gen);
      Y = GenerateY(gen);
      Z = GenerateZ(gen);
      if ((X > minX_user) and (X < maxX_user) and (Y > minY_user) and (Y < maxY_user) and (Z > minZ_user) and (Z < maxZ_user))
      {
        NonZeroCases++;
      }
    }
    MonteCarloProbability = NonZeroCases/NumSample;

    return;
  }

private:
  vtkm::Vec<vtkm::Float64, 3> InputBottomLeft;
  vtkm::Vec<vtkm::Float64, 3> InputTopRight;
};

}
}
}

#endif // vtk_m_worklet_uncertainty_Fiber_h
