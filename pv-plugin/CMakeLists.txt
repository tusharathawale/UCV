paraview_add_plugin(UncertainContour
  VERSION "0.0"
  SERVER_MANAGER_XML UncertainContour.xml
  MODULES UncertainContourVTKFilters
  MODULE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/UncertainContour/vtk.module
  )
