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