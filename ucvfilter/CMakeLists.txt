set(headers
  CriticalPointsHistogram.h
  EnsembleFieldToHistogramField.h
  )

set(sources
  CriticalPointsHistogram.cxx
  EnsembleFieldToHistogramField.cxx
  )

add_library(vtkm_filter_uncertainty
  ${sources}
  ${headers}
  )
vtkm_add_target_information(vtkm_filter_uncertainty
  DEVICE_SOURCES ${sources}
  DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
  )
target_include_directories(vtkm_filter_uncertainty
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  )
target_link_libraries(vtkm_filter_uncertainty
  PUBLIC vtkm::filter_core
  )
