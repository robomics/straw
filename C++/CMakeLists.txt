cmake_minimum_required(VERSION 3.14) # CMake version check
project(strawC) # Create project "simple_example"
set(CMAKE_CXX_STANDARD 14) # Enable c++14 standard

option(BUILD_EXAMPLES "Build straw example(s)" OFF)
add_library(straw_project_options INTERFACE)
target_compile_features(straw_project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_subdirectory(src)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
