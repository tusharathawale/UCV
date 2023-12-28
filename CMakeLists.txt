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

#Fine eigen assume this is installed by apt-get install
#find_package (Eigen3 3.3 REQUIRED)

OPTION (USE_CUDA OFF)
OPTION (USE_HIP OFF)
message ("USE_CUDA is setting as: " ${USE_CUDA})
message ("USE_HIP is setting as: " ${USE_HIP})

add_executable(LoadEnsAndProcess LoadEnsAndProcess.cpp)
target_link_libraries(LoadEnsAndProcess ${VTKm_LIBRARIES})



