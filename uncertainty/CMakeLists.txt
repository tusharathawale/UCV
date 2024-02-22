##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

set(uncertainty_headers
  Fiber.h
  )
set(uncertainty_sources_device
  Fiber.cxx
  )

vtkm_library(
  NAME vtkm_filter_uncertainty
  HEADERS ${uncertainty_headers}
  DEVICE_SOURCES ${uncertainty_sources_device}
  USE_VTKM_JOB_POOL
)

set_property(TARGET
   vtkm_filter_uncertainty
   PROPERTY MODIFY_CUDA_FLAGS TRUE
   )

vtkm_add_target_information(vtkm_filter_uncertainty
   DEVICE_SOURCES ContourUncertainUniform.cxx
   DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
  )

# vtkm_add_instantiations(instantiation_files
#   INSTANTIATIONS_FILE worklet/EntropyUniform.h
#   INSTANTIATIONS_FILE ContourUncertainUniform.h
#   )
# list(APPEND uncertainty_sources_device ${instantiation_files})

target_link_libraries(vtkm_filter PUBLIC INTERFACE vtkm_filter_uncertainty)

add_subdirectory(worklet)
