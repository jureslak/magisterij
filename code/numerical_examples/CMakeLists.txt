cmake_minimum_required(VERSION 2.8.12)
project(convergence_test)

include_directories("${CMAKE_SOURCE_DIR}/../e62numcodes/src")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
SET(CMAKE_CXX_FLAGS "-O3 -std=c++14 -DNDEBUG -DEIGEN_NO_DEBUG")
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

add_subdirectory(${CMAKE_SOURCE_DIR}/../e62numcodes/ build_util)

# be lazy and use all utils libs, you can specify one by one
add_executable(convergence_test_1D_neumann convergence_test_1D_neumann.cpp)
target_link_libraries(convergence_test_1D_neumann pde_utils)

add_executable(convergence_test_1D_rbf convergence_test_1D_rbf.cpp)
target_link_libraries(convergence_test_1D_rbf pde_utils)

add_executable(poisson_square_implicit poisson_square_implicit.cpp)
target_link_libraries(poisson_square_implicit pde_utils)

add_executable(poisson_square_implicit_basis poisson_square_implicit_basis.cpp)
target_link_libraries(poisson_square_implicit_basis pde_utils)

add_executable(poisson_square_implicit_sigma_scan poisson_square_implicit_sigma_scan.cpp)
target_link_libraries(poisson_square_implicit_sigma_scan pde_utils)

add_executable(poisson_square_explicit poisson_square_explicit.cpp)
target_link_libraries(poisson_square_explicit pde_utils)

add_executable(poisson_weird poisson_weird.cpp)
target_link_libraries(poisson_weird pde_utils)

add_executable(poisson_neumann poisson_neumann.cpp)
target_link_libraries(poisson_neumann pde_utils)
