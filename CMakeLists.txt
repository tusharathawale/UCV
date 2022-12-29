cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::rendering vtkm::source vtkm::vtkmdiympi vtkm::io)

#Find MPI
find_package(MPI REQUIRED)

add_executable(ucv_extract ucv_extract.cpp)
target_link_libraries(ucv_extract ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_umc ucv_umc.cpp)
target_link_libraries(ucv_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

#add_executable(ucv_extract_vector ucv_extract_vector.cpp)
#target_link_libraries(ucv_extract_vector ${VTKm_LIBRARIES} MPI::MPI_CXX)

#add_executable(ucv_umc_vector ucv_umc_vector.cpp)
#target_link_libraries(ucv_umc_vector ${VTKm_LIBRARIES} MPI::MPI_CXX)

#add_executable(point_to_cell point_to_cell.cxx)
#target_link_libraries(point_to_cell ${VTKm_LIBRARIES})


