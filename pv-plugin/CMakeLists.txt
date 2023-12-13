paraview_add_plugin(UncertainAnalysis
  VERSION "0.0"
  SERVER_MANAGER_XML UncertainAnalysis.xml
  MODULES UncertainAnalysisVTKFilters
  MODULE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/UncertainAnalysis/vtk.module
  )