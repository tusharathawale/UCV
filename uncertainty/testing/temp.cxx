#include <vtkm/cont/Initialize.h>
#include <vtkm/filter/contour/Slice.h>
#include <vtkm/io/VTKDataSetReader.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/source/Wavelet.h>

int main(int argc, char** argv)
{
  vtkm::cont::Initialize(argc, argv);

  vtkm::source::Wavelet wavelet;
  wavelet.SetExtent(vtkm::Id3(-8), vtkm::Id3(8));
  auto ds = wavelet.Execute();

  vtkm::Plane plane(vtkm::Plane::Vector{ 1, 1, 1 });
  vtkm::filter::contour::Slice slice;
  slice.SetImplicitFunction(plane);
  auto result = slice.Execute(ds);

  result.PrintSummary(std::cout);

  return 0;
}
