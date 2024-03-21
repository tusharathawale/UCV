paraview_add_plugin(UncertainCriticalPoints
  VERSION "0.0"
  SERVER_MANAGER_XML UncertainCriticalPoints.xml
  MODULES UncertainCriticalPointsVTKFilters
  MODULE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/UncertainCriticalPoints/vtk.module
  )
