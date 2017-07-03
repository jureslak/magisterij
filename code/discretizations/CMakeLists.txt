cmake_minimum_required(VERSION 2.8.12)
project(discretization_test)
include_directories(${CMAKE_SOURCE_DIR}/../e62numcodes/src/)
add_subdirectory(${CMAKE_SOURCE_DIR}/../e62numcodes/ build_util)

SET(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -fopenmp")

# be lazy and use all utils libs, you can specify one by one
add_executable(circle discretization_circle.cpp)
target_link_libraries(circle pde_utils)

add_executable(relax discretization_relax_example.cpp)
target_link_libraries(relax pde_utils X11)
