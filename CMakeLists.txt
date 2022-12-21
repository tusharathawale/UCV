cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED)
set(VTKm_LIBRARIES vtkm_cont vtkm_filter vtkm_rendering vtkm_source vtkmdiympi)

#Find MPI
find_package(MPI REQUIRED)

add_executable(ucv ucv.cpp)
target_link_libraries(ucv ${VTKm_LIBRARIES} MPI::MPI_CXX)
