cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES C CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED)
#set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)
#NOMPI version
message("debug vtkm_source:" ${vtkm_source})
message("debug VTKm_DIR ${VTKm_DIR}")
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::io)


#SET(CMAKE_CXX_FLAGS "-DUSE_LOG")


#Find MPI
find_package(MPI REQUIRED)

OPTION (USE_HIP OFF)

if(USE_HIP)

message ("Build messages: Using HIP GPU to build the project%%%")

set_source_files_properties(LoadEnsAndProcessByUniform.cpp PROPERTIES LANGUAGE "HIP")
add_executable(LoadEnsAndProcessByUniform LoadEnsAndProcessByUniform.cpp)
target_link_libraries(LoadEnsAndProcessByUniform ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(LoadEnsAndProcessByHistogram.cpp PROPERTIES LANGUAGE "HIP")
add_executable(LoadEnsAndProcessByHistogram LoadEnsAndProcessByHistogram.cpp)
target_link_libraries(LoadEnsAndProcessByHistogram ${VTKm_LIBRARIES} MPI::MPI_CXX)

elseif(USE_CUDA)

set_source_files_properties(./ucvworklet/CriticalPointHistogramWorklet.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(./ucvworklet/CriticalPointWorklet.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(./ucvworklet/CriticalPointWorklet.hpp PROPERTIES LANGUAGE "CUDA")


add_library(LoadEnsAndProcessByHistogram LoadEnsAndProcessByHistogram.cpp)
target_link_libraries(LoadEnsAndProcessByHistogram ${VTKm_LIBRARIES})
set_target_properties(LoadEnsAndProcessByHistogram PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

else()

add_executable(LoadEnsAndProcessByUniform LoadEnsAndProcessByUniform.cpp)
target_link_libraries(LoadEnsAndProcessByUniform ${VTKm_LIBRARIES})

add_executable(LoadEnsAndProcessByHistogram LoadEnsAndProcessByHistogram.cpp)
target_link_libraries(LoadEnsAndProcessByHistogram ${VTKm_LIBRARIES})

add_executable(LoadWeatherProcessByUniform LoadWeatherProcessByUniform.cpp)
target_link_libraries(LoadWeatherProcessByUniform ${VTKm_LIBRARIES})




endif()

