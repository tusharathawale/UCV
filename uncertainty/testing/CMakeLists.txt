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

else()

add_executable(TestFiber TestFiber.cxx ../Fiber.cxx)
target_link_libraries(TestFiber ${VTKm_LIBRARIES})

add_executable(TestSuperNova TestSuperNova.cxx ../Fiber.cxx)
target_link_libraries(TestSuperNova ${VTKm_LIBRARIES})

add_executable(FiberSurfaceExample FiberSurfaceExample.cxx ../Fiber.cxx)
target_link_libraries(FiberSurfaceExample ${VTKm_LIBRARIES})

add_executable(DecomposeSuperNovaToEns DecomposeSuperNovaToEns.cxx)
target_link_libraries(DecomposeSuperNovaToEns ${VTKm_LIBRARIES})

endif()
