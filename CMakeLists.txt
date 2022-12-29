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

add_executable(ucv_extract_vector ucv_extract_vector.cpp)
target_link_libraries(ucv_extract_vector ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_umc_vector ucv_umc_vector.cpp)
target_link_libraries(ucv_umc_vector ${VTKm_LIBRARIES} MPI::MPI_CXX)

#Not sure why it takes so long to build this
#add_executable(ucv_umc ucv_umc.cpp)
#target_link_libraries(ucv_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(point_to_cell point_to_cell.cxx)
target_link_libraries(point_to_cell ${VTKm_LIBRARIES})

#add_executable(point_to_cell point_to_cell.cxx)
#target_link_libraries(point_to_cell vtkm_cont vtkm_filter vtkm_io)
#vtkm_add_target_information(point_to_cell
#  DROP_UNUSED_SYMBOLS MODIFY_CUDA_FLAGS
#  DEVICE_SOURCES point_to_cell.cxx)

