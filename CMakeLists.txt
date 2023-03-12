cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)
set(VTKm_LIBRARIESNOMPI vtkm::cont vtkm::filter vtkm::source vtkm::io)

#Find MPI
find_package(MPI REQUIRED)

#Fine eigen assume this is installed by apt-get install
#find_package (Eigen3 3.3 REQUIRED)

OPTION (USE_CUDA OFF)
OPTION (USE_HIP OFF)
message ("USE_CUDA is setting as: " ${USE_CUDA})
message ("USE_HIP is setting as: " ${USE_HIP})

if(USE_CUDA)
set_source_files_properties(ucv_reduce_umc.cpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/EntropyUniform.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanRaw.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanStdev.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMinMax.hpp PROPERTIES LANGUAGE "CUDA")

add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
set_target_properties(ucv_reduce_umc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

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

set_source_files_properties(test_mvg_supernova.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvg_supernova test_mvg_supernova.cpp)
set_target_properties(test_mvg_supernova PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvg_supernova ${VTKm_LIBRARIES})

set_source_files_properties(uvm_point_neighborhood.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(uvm_point_neighborhood uvm_point_neighborhood.cpp)
set_target_properties(uvm_point_neighborhood PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(uvm_point_neighborhood ${VTKm_LIBRARIES})

elseif(USE_HIP)

set_source_files_properties(ucvworklet/EntropyUniform.hpp PROPERTIES LANGUAGE "HIP")
set_source_files_properties(ucvworklet/ExtractingMeanRaw.hpp PROPERTIES LANGUAGE "HIP")
set_source_files_properties(ucvworklet/ExtractingMeanStdev.hpp PROPERTIES LANGUAGE "HIP")
set_source_files_properties(ucvworklet/ExtractingMinMax.hpp PROPERTIES LANGUAGE "HIP")

set_source_files_properties(ucv_reduce_umc.cpp PROPERTIES LANGUAGE "HIP")
add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIESNOMPI} MPI::MPI_CXX)

set_source_files_properties(test_mvgaussian_wind.cpp PROPERTIES LANGUAGE "HIP")
add_executable(test_mvgaussian_wind test_mvgaussian_wind.cpp)
target_link_libraries(test_mvgaussian_wind ${VTKm_LIBRARIESNOMPI})

set_source_files_properties(test_mvgaussian_redsea.cpp PROPERTIES LANGUAGE "HIP")
add_executable(test_mvgaussian_redsea test_mvgaussian_redsea.cpp)
target_link_libraries(test_mvgaussian_redsea ${VTKm_LIBRARIESNOMPI})

set_source_files_properties(uvm_point_neighborhood.cpp PROPERTIES LANGUAGE "HIP")
add_executable(uvm_point_neighborhood uvm_point_neighborhood.cpp)
target_link_libraries(uvm_point_neighborhood ${VTKm_LIBRARIESNOMPI})

else()
 
add_executable(ucv_extract ucv_extract.cpp)
target_link_libraries(ucv_extract ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_extract_reduce ucv_extract_reduce.cpp)
target_link_libraries(ucv_extract_reduce ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_umc ucv_umc.cpp)
target_link_libraries(ucv_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

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

add_executable(test_mvg_supernova test_mvg_supernova.cpp)
target_link_libraries(test_mvg_supernova ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(compute_cases compute_cases.cpp)
target_link_libraries(compute_cases)

add_executable(test_ucv_matrix ./ucvworklet/linalg/test_ucv_matrix.cpp)
target_link_libraries(test_ucv_matrix ${VTKm_LIBRARIES})

add_executable(uvm_point_neighborhood uvm_point_neighborhood.cpp)
target_link_libraries(uvm_point_neighborhood ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_3by3 ./ucvworklet/linalg/test_ucv_matrix_static_3by3.cpp)
target_link_libraries(test_ucv_matrix_static_3by3 ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_4by4 ./ucvworklet/linalg/test_ucv_matrix_static_4by4.cpp)
target_link_libraries(test_ucv_matrix_static_4by4 ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_8by8 ./ucvworklet/linalg/test_ucv_matrix_static_8by8.cpp)
target_link_libraries(test_ucv_matrix_static_8by8 ${VTKm_LIBRARIES})

endif()



