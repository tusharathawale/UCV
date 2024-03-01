paraview_add_plugin(UncertainFiberSurface
  VERSION "0.0"
  SERVER_MANAGER_XML UncertainFiberSurface.xml
  MODULES UncertainFiberSurfaceVTKFilters
  MODULE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/UncertainFiberSurface/vtk.module
  )
