set(classes
  vtkCriticalPointsHistogram
  vtkEnsembleFieldToHistogramField
  )

vtk_module_add_module(UncertainCriticalPointsVTKFilters
  CLASSES ${classes}
  )

target_sources(UncertainCriticalPointsVTKFilters
  PRIVATE
    VTKToVTKm.h
    VTKToVTKm.cxx
)
vtk_module_link(UncertainCriticalPointsVTKFilters
  PRIVATE vtkm_filter_uncertainty
  )
vtkm_add_target_information(vtkm_filter_uncertainty
  DEVICE_SOURCES VTKToVTKm.cxx
  DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
  )
