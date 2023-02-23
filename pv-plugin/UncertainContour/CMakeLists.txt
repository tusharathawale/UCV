set(classes
  vtkContourUncertainEnsemble
  vtkContourUncertainIndependentGaussian
  vtkContourUncertainUniform
  vtkSubsampleUncertaintyEnsemble
  vtkSubsampleUncertaintyIndependentGaussian
  vtkSubsampleUncertaintyUniform
  )

vtk_module_add_module(UncertainContourVTKFilters
  CLASSES ${classes}
  )
vtk_module_link(UncertainContourVTKFilters
  PRIVATE filter_uncertainty
  )
