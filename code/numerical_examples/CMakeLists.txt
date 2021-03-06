cmake_minimum_required(VERSION 2.8.12)

#SET(CMAKE_C_COMPILER "/opt/intel/compilers_and_libraries_2016.1.150/linux/bin/intel64/icc")
#SET(CMAKE_CXX_COMPILER "/opt/intel/compilers_and_libraries_2016.1.150/linux/bin/intel64/icpc")

project(numerical_examples)

include_directories("${CMAKE_SOURCE_DIR}/../e62numcodes/src")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
SET(CMAKE_CXX_FLAGS "-O3 -std=c++11")
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

add_subdirectory(${CMAKE_SOURCE_DIR}/../e62numcodes/ build_util)

# be lazy and use all utils libs, you can specify one by one
add_executable(convergence_test_1D_neumann convergence_test_1D_neumann.cpp)
target_compile_definitions(convergence_test_1D_neumann PUBLIC NDEBUG)
target_link_libraries(convergence_test_1D_neumann pde_utils)

add_executable(convergence_test_1D_rbf convergence_test_1D_rbf.cpp)
target_link_libraries(convergence_test_1D_rbf pde_utils)

add_executable(poisson_square_implicit poisson_square_implicit.cpp)
target_compile_definitions(poisson_square_implicit PUBLIC NDEBUG)
target_link_libraries(poisson_square_implicit pde_utils)

add_executable(poisson_square_implicit_rbf_konv_k_mon poisson_square_implicit_rbf_konv_k_mon.cpp)
target_link_libraries(poisson_square_implicit_rbf_konv_k_mon pde_utils)

add_executable(poisson_square_implicit_sigma_scan poisson_square_implicit_sigma_scan.cpp)
target_compile_definitions(poisson_square_implicit_sigma_scan PUBLIC NDEBUG)
target_link_libraries(poisson_square_implicit_sigma_scan pde_utils)

add_executable(poisson_square_explicit poisson_square_explicit.cpp)
target_link_libraries(poisson_square_explicit pde_utils)

#add_executable(poisson_square_parallel poisson_square_parallel.cpp)
#target_compile_options(poisson_square_parallel PUBLIC "-openmp")
#target_compile_definitions(poisson_square_parallel PUBLIC EIGEN_USE_MKL_ALL MKL_LP64 NDEBUG)
#link_directories("$ENV{MKLROOT}/lib/intel64")
#target_link_libraries(poisson_square_parallel pde_utils mkl_intel_lp64 mkl_intel_thread mkl_core pthread iomp5)

add_executable(poisson_weird poisson_weird.cpp)
target_link_libraries(poisson_weird pde_utils)

add_executable(poisson_neumann poisson_neumann.cpp)
target_link_libraries(poisson_neumann pde_utils)

add_executable(hertzian_convergence hertzian_convergence.cpp overseer_hertzian.hpp)
target_compile_definitions(hertzian_convergence PUBLIC NDEBUG)
target_link_libraries(hertzian_convergence pde_utils)

add_executable(hertzian_refined_convergence hertzian_refine_levels_convergence.cpp overseer_hertzian.hpp)
target_compile_definitions(hertzian_refined_convergence PUBLIC NDEBUG)
target_link_libraries(hertzian_refined_convergence pde_utils)

add_executable(hertzian_matrix_example hertzian_matrix_example.cpp overseer_hertzian.hpp)
target_link_libraries(hertzian_matrix_example pde_utils)

add_executable(hertzian_domain_too_small hertzian_domain_too_small.cpp overseer_hertzian.hpp)
target_link_libraries(hertzian_domain_too_small pde_utils)

add_executable(cantilever_beam cantilever_beam.cpp overseer_cantilever.hpp)
target_compile_definitions(cantilever_beam PUBLIC NDEBUG)
target_link_libraries(cantilever_beam pde_utils)

add_executable(cantilever_beam_matrix_example cantilever_beam_matrix_example.cpp overseer_cantilever.hpp)
target_link_libraries(cantilever_beam_matrix_example pde_utils)

add_executable(realistic_beam realistic_beam.cpp overseer_cantilever.hpp)
target_link_libraries(realistic_beam pde_utils)

add_executable(cantilever_beam_with_holes cantilever_beam_with_holes.cpp overseer_cantilever.hpp)
target_link_libraries(cantilever_beam_with_holes pde_utils)

add_executable(fwo_convergence fwo_convergence.cpp overseer_fwo.hpp)
target_link_libraries(fwo_convergence pde_utils)

add_executable(fwo_cases fwo_cases.cpp overseer_fwo.hpp)
target_link_libraries(fwo_cases pde_utils)

add_executable(fwo_table fwo_table.cpp overseer_fwo.hpp)
target_link_libraries(fwo_table pde_utils)

add_executable(fwo_solution fwo_solution.cpp overseer_fwo.hpp)
target_link_libraries(fwo_solution pde_utils)
