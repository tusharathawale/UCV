cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES C CXX)
add_compile_options(-std=c++17)

option(BUILD_PARAVIEW_PLUGIN
  "Build a plugin for ParaView that makes uncertainty filters available."
  OFF
  )
if(BUILD_PARAVIEW_PLUGIN)
  find_package(ParaView REQUIRED)
endif()

#Find the VTK-m package
find_package(VTKm REQUIRED)
#set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)
#NOMPI version
message("debug vtkm_source:" ${vtkm_source})
message("debug VTKm_DIR ${VTKm_DIR}")
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::io vtkm_filter_uncertainty)


#SET(CMAKE_CXX_FLAGS "-DUSE_LOG")


#Find MPI

add_executable(LoadEnsAndProcessByHistogram LoadEnsAndProcessByHistogram.cpp)
target_link_libraries(LoadEnsAndProcessByHistogram ${VTKm_LIBRARIES})

OPTION (USE_HIP OFF)

add_subdirectory(ucvfilter)

if(USE_HIP)

message ("Build messages: Using HIP GPU to build the project%%%")

set_source_files_properties(LoadEnsAndProcessByUniform.cpp PROPERTIES LANGUAGE "HIP")
add_executable(LoadEnsAndProcessByUniform LoadEnsAndProcessByUniform.cpp)
target_link_libraries(LoadEnsAndProcessByUniform ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(LoadWeatherProcessByUniform.cpp PROPERTIES LANGUAGE "HIP")
add_executable(LoadWeatherProcessByUniform LoadWeatherProcessByUniform.cpp)
target_link_libraries(LoadWeatherProcessByUniform ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(LoadWeatherProcessByMC.cpp PROPERTIES LANGUAGE "HIP")
add_executable(LoadWeatherProcessByMC LoadWeatherProcessByMC.cpp)
target_link_libraries(LoadWeatherProcessByMC ${VTKm_LIBRARIES} MPI::MPI_CXX)

elseif(USE_CUDA)

set_source_files_properties(./ucvworklet/CriticalPointMonteCarlo.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(./ucvworklet/CriticalPointWorklet.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(./ucvworklet/CriticalPointWorkletAvoidOverflow.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(./ucvworklet/ExtractMinMaxByErr.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(./ucvworklet/ExtractMinMaxOfPoint.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(LoadWeatherProcessByMC.cpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(LoadWeatherProcessByUniform.cpp PROPERTIES LANGUAGE "CUDA")

add_executable(LoadWeatherProcessByMC LoadWeatherProcessByMC.cpp)
target_link_libraries(LoadWeatherProcessByMC ${VTKm_LIBRARIES})
set_target_properties(LoadWeatherProcessByMC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(LoadWeatherProcessByUniform LoadWeatherProcessByUniform.cpp)
target_link_libraries(LoadWeatherProcessByUniform ${VTKm_LIBRARIES})
set_target_properties(LoadWeatherProcessByUniform PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

else()

add_executable(LoadEnsAndProcessByUniform LoadEnsAndProcessByUniform.cpp)
target_link_libraries(LoadEnsAndProcessByUniform ${VTKm_LIBRARIES})

add_executable(LoadWeatherProcessByUniform LoadWeatherProcessByUniform.cpp)
target_link_libraries(LoadWeatherProcessByUniform ${VTKm_LIBRARIES})

add_executable(LoadEnsAndProcessByMC LoadEnsAndProcessByMC.cpp)
target_link_libraries(LoadEnsAndProcessByMC ${VTKm_LIBRARIES})

add_executable(LoadWeatherProcessByMC LoadWeatherProcessByMC.cpp)
target_link_libraries(LoadWeatherProcessByMC ${VTKm_LIBRARIES})



endif()

if(BUILD_PARAVIEW_PLUGIN)
  # ParaView plugins need this included to be built.
  # (Why don't they just include themselves?)
  include(GNUInstallDirs)

  # ParaView plugins need to be packaged in a shared library
  set(BUILD_SHARED_LIBS ON)

  paraview_plugin_scan(
    PLUGIN_FILES ${CMAKE_CURRENT_SOURCE_DIR}/pv-plugin/paraview.plugin
    PROVIDES_PLUGINS plugins
    ENABLE_BY_DEFAULT ON
    HIDE_PLUGINS_FROM_CACHE ON
    )

  # Note: this will call subdirectory for pv-plugin so we don't have to.
  paraview_plugin_build(PLUGINS ${plugins})
endif()
