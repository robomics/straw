cmake_minimum_required(VERSION 3.13)
cmake_policy(VERSION 3.13...3.27)

project(strawC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)


add_library(straw_project_options INTERFACE)
target_compile_features(straw_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

add_subdirectory(src)
