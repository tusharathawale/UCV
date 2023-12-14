set(classes
    vtkUncertainContour
  )

vtk_module_add_module(UncertainAnalysisVTKFilters
  CLASSES ${classes}
  )

#if the inner code is built by separate filter
#using this pattern
#vtk_module_link(UncertainAnalysisVTKFilters
#  PRIVATE filter_uncertainty
# )

find_package(VTK COMPONENTS
IOLegacy
IOXML
AcceleratorsVTKmCore
AcceleratorsVTKmDataModel
CommonExecutionModel
)

add_executable(TestToVTKm TestToVTKm.cpp)
target_link_libraries(TestToVTKm ${VTKm_LIBRARIES} ${VTK_LIBRARIES})