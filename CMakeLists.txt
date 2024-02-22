cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::io)

add_subdirectory(uncertainty)