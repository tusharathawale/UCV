cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::rendering vtkm::source vtkm::vtkmdiympi vtkm::io)

#Find MPI
find_package(MPI REQUIRED)

#Fine eigen assume this is installed by apt-get install
find_package (Eigen3 3.3 REQUIRED)

OPTION (USE_GPU OFF)
message ("USE_GPU is setting as: " ${USE_GPU})
if(USE_GPU)

set_source_files_properties(ucv_reduce_umc.cpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/EntropyUniform.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanRaw.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanStdev.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMinMax.hpp PROPERTIES LANGUAGE "CUDA")

add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
set_target_properties(ucv_reduce_umc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX Eigen3::Eigen)

set_source_files_properties(test_mvgaussian_wind.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvgaussian_wind test_mvgaussian_wind.cpp)
set_target_properties(test_mvgaussian_wind PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvgaussian_wind ${VTKm_LIBRARIES})


add_executable(test_ucv_matrix ./ucvworklet/linalg/test_ucv_matrix.cpp)
target_link_libraries(test_ucv_matrix ${VTKm_LIBRARIES})

else()
 
add_executable(ucv_extract ucv_extract.cpp)
target_link_libraries(ucv_extract ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_extract_reduce ucv_extract_reduce.cpp)
target_link_libraries(ucv_extract_reduce ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_umc ucv_umc.cpp)
target_link_libraries(ucv_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX Eigen3::Eigen)

add_executable(test_mvgaussian_wind test_mvgaussian_wind.cpp)
target_link_libraries(test_mvgaussian_wind ${VTKm_LIBRARIES} MPI::MPI_CXX Eigen3::Eigen)

add_executable(test_mvgaussian_3d test_mvgaussian_3d.cpp)
target_link_libraries(test_mvgaussian_3d ${VTKm_LIBRARIES} MPI::MPI_CXX Eigen3::Eigen)

add_executable(compute_cases compute_cases.cpp)
target_link_libraries(compute_cases)

add_executable(test_ucv_matrix ./ucvworklet/linalg/test_ucv_matrix.cpp)
target_link_libraries(test_ucv_matrix ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_4by4 ./ucvworklet/linalg/test_ucv_matrix_static_4by4.cpp)
target_link_libraries(test_ucv_matrix_static_4by4 ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_8by8 ./ucvworklet/linalg/test_ucv_matrix_static_8by8.cpp)
target_link_libraries(test_ucv_matrix_static_8by8 ${VTKm_LIBRARIES})

endif()



