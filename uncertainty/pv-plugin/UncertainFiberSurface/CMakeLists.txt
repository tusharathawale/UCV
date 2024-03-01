set(classes
  vtkFiberSurfaceUncertain
  )

vtk_module_add_module(UncertainFiberSurfaceVTKFilters
  CLASSES ${classes}
  )
vtk_module_link(UncertainFiberSurfaceVTKFilters
  PRIVATE filter_uncertainty
  )
