cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED)
#set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)
#NOMPI version
message("debug vtkm_source:" ${vtkm_source})
message("debug VTKm_DIR ${VTKm_DIR}")
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::io)


#Find MPI
find_package(MPI REQUIRED)
add_executable(LoadEnsAndProcessByUniform LoadEnsAndProcessByUniform.cpp)
target_link_libraries(LoadEnsAndProcessByUniform ${VTKm_LIBRARIES})

add_executable(LoadEnsAndProcessByHistogram LoadEnsAndProcessByHistogram.cpp)
target_link_libraries(LoadEnsAndProcessByHistogram ${VTKm_LIBRARIES})



