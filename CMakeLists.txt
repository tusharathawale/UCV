cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX C)
add_compile_options(-std=c++14)

option(BUILD_PARAVIEW_PLUGIN
  "Build a plugin for ParaView that makes uncertainty filters available."
  OFF
  )
if(BUILD_PARAVIEW_PLUGIN)
  find_package(ParaView REQUIRED)
endif()

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)
message("${VTKm_DIR}")

#Find MPI
find_package(MPI REQUIRED)

#Fine eigen assume this is installed by apt-get install
#find_package (Eigen3 3.3 REQUIRED)

set(filter_sources
  ContourUncertainEnsemble.cxx
  ContourUncertainIndependentGaussian.cxx
  ContourUncertainUniform.cxx
  SubsampleUncertaintyEnsemble.cxx
  SubsampleUncertaintyIndependentGaussian.cxx
  SubsampleUncertaintyUniform.cxx
  )

OPTION (USE_GPU "Compile GPU support." OFF)
message ("USE_GPU is setting as: " ${USE_GPU})
if(USE_GPU)

set_source_files_properties(ucv_reduce_umc.cpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/EntropyUniform.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanRaw.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanStdev.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMinMax.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(${filter_sources} PROPERTIES LANGAUGE "CUDA")

add_library(filter_uncertainty ${filter_sources})
target_link_libraries(filter_uncertainty PUBLIC ${VTKm_LIBRARIES})

add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
set_target_properties(ucv_reduce_umc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX filter_uncertainty)

set_source_files_properties(test_mvgaussian_wind.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvgaussian_wind test_mvgaussian_wind.cpp)
set_target_properties(test_mvgaussian_wind PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvgaussian_wind ${VTKm_LIBRARIES})

set_source_files_properties(test_mvgaussian_redsea_mpi.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvgaussian_redsea_mpi test_mvgaussian_redsea_mpi.cpp)
set_target_properties(test_mvgaussian_redsea_mpi PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvgaussian_redsea_mpi ${VTKm_LIBRARIES})

set_source_files_properties(test_mvgaussian_redsea.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvgaussian_redsea test_mvgaussian_redsea.cpp)
set_target_properties(test_mvgaussian_redsea PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvgaussian_redsea ${VTKm_LIBRARIES})

#add_executable(test_ucv_matrix ./ucvworklet/linalg/test_ucv_matrix.cpp)
#target_link_libraries(test_ucv_matrix ${VTKm_LIBRARIES})

set_source_files_properties(test_mvgaussian_redsea_as3d.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvgaussian_redsea_as3d test_mvgaussian_redsea_as3d.cpp)
set_target_properties(test_mvgaussian_redsea_as3d PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvgaussian_redsea_as3d ${VTKm_LIBRARIES})


else()

add_library(filter_uncertainty ${filter_sources})
target_link_libraries(filter_uncertainty PUBLIC ${VTKm_LIBRARIES})
target_include_directories(filter_uncertainty
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  )

add_executable(ucv_extract ucv_extract.cpp)
target_link_libraries(ucv_extract ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_extract_reduce ucv_extract_reduce.cpp)
target_link_libraries(ucv_extract_reduce ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_umc ucv_umc.cpp)
target_link_libraries(ucv_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX filter_uncertainty)

add_executable(test_mvgaussian_wind test_mvgaussian_wind.cpp)
target_link_libraries(test_mvgaussian_wind ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_redsea test_mvgaussian_redsea.cpp)
target_link_libraries(test_mvgaussian_redsea ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_redsea_mpi test_mvgaussian_redsea_mpi.cpp)
target_link_libraries(test_mvgaussian_redsea_mpi ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_redsea_as3d test_mvgaussian_redsea_as3d.cpp)
target_link_libraries(test_mvgaussian_redsea_as3d ${VTKm_LIBRARIES} MPI::MPI_CXX)

#add_executable(test_mvgaussian_3d test_mvgaussian_3d.cpp)
#target_link_libraries(test_mvgaussian_3d ${VTKm_LIBRARIES} MPI::MPI_CXX Eigen3::Eigen)

add_executable(compute_cases compute_cases.cpp)
target_link_libraries(compute_cases)

add_executable(test_ucv_matrix ./ucvworklet/linalg/test_ucv_matrix.cpp)
target_link_libraries(test_ucv_matrix ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_3by3 ./ucvworklet/linalg/test_ucv_matrix_static_3by3.cpp)
target_link_libraries(test_ucv_matrix_static_3by3 ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_4by4 ./ucvworklet/linalg/test_ucv_matrix_static_4by4.cpp)
target_link_libraries(test_ucv_matrix_static_4by4 ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_8by8 ./ucvworklet/linalg/test_ucv_matrix_static_8by8.cpp)
target_link_libraries(test_ucv_matrix_static_8by8 ${VTKm_LIBRARIES})

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
