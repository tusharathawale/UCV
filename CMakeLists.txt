cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

option(BUILD_PARAVIEW_PLUGIN
  "Build a plugin for ParaView that makes uncertainty filters available."
  OFF
  )
if(BUILD_PARAVIEW_PLUGIN)
  find_package(ParaView REQUIRED)
endif()

#Find the VTK-m package
find_package(VTKm REQUIRED)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::io)

add_subdirectory(uncertainty)
