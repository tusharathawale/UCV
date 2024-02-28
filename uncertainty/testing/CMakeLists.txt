##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

OPTION (USE_HIP OFF)

if(USE_HIP)

set_source_files_properties(TestFiber.cxx PROPERTIES LANGUAGE "HIP")
add_executable(TestFiber TestFiber.cxx ../Fiber.cxx)
target_link_libraries(TestFiber ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(TestSuperNova.cxx PROPERTIES LANGUAGE "HIP")
add_executable(TestSuperNova TestSuperNova.cxx ../Fiber.cxx)
target_link_libraries(TestSuperNova ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(TestSuperNovaComparison.cxx PROPERTIES LANGUAGE "HIP")
add_executable(TestSuperNovaComparison TestSuperNovaComparison.cxx ../Fiber.cxx)
target_link_libraries(TestSuperNovaComparison ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(TestRedSea.cxx PROPERTIES LANGUAGE "HIP")
add_executable(TestRedSea TestRedSea.cxx ../Fiber.cxx)
target_link_libraries(TestRedSea ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(TestRedSeaComparison.cxx PROPERTIES LANGUAGE "HIP")
add_executable(TestRedSeaComparison TestRedSeaComparison.cxx ../Fiber.cxx)
target_link_libraries(TestRedSeaComparison ${VTKm_LIBRARIES} MPI::MPI_CXX)

elseif(USE_CUDA)

message("---use cuda build---")

set_source_files_properties(../Fiber.h PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(../worklet/Fiber.h PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(../worklet/ExtractingByNeigoborhoodEnsTwoFields.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(../worklet/ExtractingMinMax.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(TestSuperNova.cxx PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(TestRedSea.cxx PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(../Fiber.cxx PROPERTIES LANGUAGE "CUDA")

add_library(fiber ../Fiber.cxx)
target_link_libraries(fiber ${VTKm_LIBRARIES})
set_target_properties(fiber PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(TestSuperNova TestSuperNova.cxx)
target_link_libraries(TestSuperNova fiber)
set_target_properties(TestSuperNova PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
# DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON can fix multiple defination issue

add_executable(TestRedSea TestRedSea.cxx)
target_link_libraries(TestRedSea fiber)
set_target_properties(TestRedSea PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

else()

add_executable(TestFiber TestFiber.cxx ../Fiber.cxx)
target_link_libraries(TestFiber ${VTKm_LIBRARIES})

add_executable(TestSuperNova TestSuperNova.cxx ../Fiber.cxx)
target_link_libraries(TestSuperNova ${VTKm_LIBRARIES})

add_executable(TestSuperNovaComparison TestSuperNovaComparison.cxx ../Fiber.cxx)
target_link_libraries(TestSuperNovaComparison ${VTKm_LIBRARIES})

add_executable(FiberSurfaceExample FiberSurfaceExample.cxx ../Fiber.cxx)
target_link_libraries(FiberSurfaceExample ${VTKm_LIBRARIES})

add_executable(DecomposeSuperNovaToEns DecomposeSuperNovaToEns.cxx)
target_link_libraries(DecomposeSuperNovaToEns ${VTKm_LIBRARIES})

add_executable(TestRedSea TestRedSea.cxx ../Fiber.cxx)
target_link_libraries(TestRedSea ${VTKm_LIBRARIES})

endif()
