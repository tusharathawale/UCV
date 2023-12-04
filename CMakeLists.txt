cmake_minimum_required(VERSION 3.12...3.15 FATAL_ERROR)
project(ucv LANGUAGES CXX)
add_compile_options(-std=c++14)

#Find the VTK-m package
find_package(VTKm REQUIRED QUIET)
#set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::vtkmdiympi vtkm::io)
#NOMPI version
set(VTKm_LIBRARIES vtkm::cont vtkm::filter vtkm::source vtkm::io)

#Find MPI
find_package(MPI REQUIRED)

#Fine eigen assume this is installed by apt-get install
#find_package (Eigen3 3.3 REQUIRED)

OPTION (USE_CUDA OFF)
OPTION (USE_HIP OFF)
message ("USE_CUDA is setting as: " ${USE_CUDA})
message ("USE_HIP is setting as: " ${USE_HIP})

if(USE_CUDA)
set_source_files_properties(ucvworklet/EntropyUniform.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanRaw.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMeanStdev.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/ExtractingMinMax.hpp PROPERTIES LANGUAGE "CUDA")
set_source_files_properties(ucvworklet/MVGaussianWithEnsemble3DTryEL.hpp PROPERTIES LANGUAGE "CUDA")

set_source_files_properties(ucv_reduce_umc.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
set_target_properties(ucv_reduce_umc PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

set_source_files_properties(ucv_reduce_umc_el.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(ucv_reduce_umc_el ucv_reduce_umc_el.cpp)
set_target_properties(ucv_reduce_umc_el PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(ucv_reduce_umc_el ${VTKm_LIBRARIES} MPI::MPI_CXX)

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

# there is compiling issues working on cuda for gsl
#find_package(GSL REQUIRED)
#set_source_files_properties(test_mvgaussian_redsea_gsl.cpp PROPERTIES LANGUAGE "CUDA")
#add_executable(test_mvgaussian_redsea_gsl test_mvgaussian_redsea_gsl.cpp)
#set_target_properties(test_mvgaussian_redsea_gsl PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#target_link_libraries(test_mvgaussian_redsea_gsl ${VTKm_LIBRARIES} GSL::gsl)

set_source_files_properties(test_mvgaussian_redsea_el.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvgaussian_redsea_el test_mvgaussian_redsea_el.cpp)
set_target_properties(test_mvgaussian_redsea_el PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvgaussian_redsea_el ${VTKm_LIBRARIES})

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

set_source_files_properties(test_mvgaussian_el_multiple_worklet.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_mvgaussian_el_multiple_worklet test_mvgaussian_el_multiple_worklet.cpp)
set_target_properties(test_mvgaussian_el_multiple_worklet PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_mvgaussian_el_multiple_worklet ${VTKm_LIBRARIES})

set_source_files_properties(test_syntheticdata_el.cpp PROPERTIES LANGUAGE "CUDA")
add_executable(test_syntheticdata_el test_syntheticdata_el.cpp)
set_target_properties(test_syntheticdata_el PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(test_syntheticdata_el ${VTKm_LIBRARIES})

elseif(USE_HIP)

set_source_files_properties(ucvworklet/EntropyUniform.hpp PROPERTIES LANGUAGE "HIP")
set_source_files_properties(ucvworklet/ExtractingMeanRaw.hpp PROPERTIES LANGUAGE "HIP")
set_source_files_properties(ucvworklet/ExtractingMeanStdev.hpp PROPERTIES LANGUAGE "HIP")
set_source_files_properties(ucvworklet/ExtractingMinMax.hpp PROPERTIES LANGUAGE "HIP")

set_source_files_properties(ucv_reduce_umc.cpp PROPERTIES LANGUAGE "HIP")
add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIESNOMPI} MPI::MPI_CXX)

set_source_files_properties(ucv_reduce_umc_el.cpp PROPERTIES LANGUAGE "HIP")
add_executable(ucv_reduce_umc_el ucv_reduce_umc_el.cpp)
target_link_libraries(ucv_reduce_umc_el ${VTKm_LIBRARIESNOMPI} MPI::MPI_CXX)

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

add_executable(ucv_extract_two_variables ucv_extract_two_variables.cpp)
target_link_libraries(ucv_extract_two_variables ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_extract_reduce ucv_extract_reduce.cpp)
target_link_libraries(ucv_extract_reduce ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_umc ucv_umc.cpp)
target_link_libraries(ucv_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_reduce_umc ucv_reduce_umc.cpp)
target_link_libraries(ucv_reduce_umc ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_reduce_umc_kde ucv_reduce_umc_kde.cpp)
target_link_libraries(ucv_reduce_umc_kde ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_reduce_umc_el ucv_reduce_umc_el.cpp)
target_link_libraries(ucv_reduce_umc_el ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(ucv_minmax_extraction ucv_minmax_extraction.cpp)
target_link_libraries(ucv_minmax_extraction ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_wind test_mvgaussian_wind.cpp)
target_link_libraries(test_mvgaussian_wind ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_redsea test_mvgaussian_redsea.cpp)
target_link_libraries(test_mvgaussian_redsea ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_redsea_gsl test_mvgaussian_redsea_gsl.cpp)
target_link_libraries(test_mvgaussian_redsea_gsl GSL::gsl ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_redsea_el test_mvgaussian_redsea_el.cpp)
target_link_libraries(test_mvgaussian_redsea_el ${VTKm_LIBRARIES} MPI::MPI_CXX)


add_executable(test_mvgaussian_redsea_mpi test_mvgaussian_redsea_mpi.cpp)
target_link_libraries(test_mvgaussian_redsea_mpi ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(test_mvgaussian_redsea_as3d test_mvgaussian_redsea_as3d.cpp)
target_link_libraries(test_mvgaussian_redsea_as3d ${VTKm_LIBRARIES} MPI::MPI_CXX)


add_executable(test_mvgaussian_redsea_entropy_diff test_mvgaussian_redsea_entropy_diff.cpp)
target_link_libraries(test_mvgaussian_redsea_entropy_diff ${VTKm_LIBRARIES} MPI::MPI_CXX)

#add_executable(test_mvgaussian_3d test_mvgaussian_3d.cpp)
#target_link_libraries(test_mvgaussian_3d ${VTKm_LIBRARIES} MPI::MPI_CXX Eigen3::Eigen)

add_executable(test_mvg_supernova test_mvg_supernova.cpp)
target_link_libraries(test_mvg_supernova ${VTKm_LIBRARIES} MPI::MPI_CXX)

add_executable(compute_cases compute_cases.cpp)
target_link_libraries(compute_cases)

add_executable(test_ucv_matrix ./ucvworklet/linalg/StaticLinalg/test_ucv_matrix.cpp)
target_link_libraries(test_ucv_matrix ${VTKm_LIBRARIES})

add_executable(uvm_point_neighborhood uvm_point_neighborhood.cpp)
target_link_libraries(uvm_point_neighborhood ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_3by3 ./ucvworklet/linalg/StaticLinalg/test_ucv_matrix_static_3by3.cpp)
target_link_libraries(test_ucv_matrix_static_3by3 ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_4by4 ./ucvworklet/linalg/StaticLinalg/test_ucv_matrix_static_4by4.cpp)
target_link_libraries(test_ucv_matrix_static_4by4 ${VTKm_LIBRARIES})

add_executable(test_ucv_matrix_static_8by8 ./ucvworklet/linalg/StaticLinalg/test_ucv_matrix_static_8by8.cpp)
target_link_libraries(test_ucv_matrix_static_8by8 ${VTKm_LIBRARIES})

find_package(GSL REQUIRED)
add_executable(sample_gsl ./ucvworklet/linalg/trygsl/sample_gsl.cpp)
target_link_libraries(sample_gsl ${VTKm_LIBRARIES} GSL::gsl)

add_executable(test_gsl_4by4 ./ucvworklet/linalg/trygsl/test_gsl_4by4.cpp)
target_link_libraries(test_gsl_4by4 GSL::gsl ${VTKm_LIBRARIES})

add_executable(test_syntheticdata_el test_syntheticdata_el.cpp)
target_link_libraries(test_syntheticdata_el GSL::gsl ${VTKm_LIBRARIES})

add_executable(test_vtkm_runvec test_vtkm_runvec.cpp)
target_link_libraries(test_vtkm_runvec ${VTKm_LIBRARIES})

add_executable(test_3D_el test_3D_el.cpp)
target_link_libraries(test_3D_el ${VTKm_LIBRARIES})

add_executable(test_syntheticdata_el_sequence_test_runtime test_syntheticdata_el_sequence_test_runtime.cpp)
target_link_libraries(test_syntheticdata_el_sequence_test_runtime ${VTKm_LIBRARIES})

add_executable(test_mvgaussian_el_multiple_worklet test_mvgaussian_el_multiple_worklet.cpp)
target_link_libraries(test_mvgaussian_el_multiple_worklet ${VTKm_LIBRARIES})

endif()



