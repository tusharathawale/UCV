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

set_source_files_properties(ucv_reduce_umc.cpp PROPERTIES LANGUAGE "HIP")
add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

else()

add_executable(TestFiber TestFiber.cxx ../Fiber.cxx)
target_link_libraries(TestFiber ${VTKm_LIBRARIES})

endif()