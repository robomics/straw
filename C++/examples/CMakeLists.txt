# Add main.cpp file of project root directory as source file
set(SOURCE_FILES examples/main.cpp src/straw.cpp)
add_executable(example example.cpp)
target_link_libraries(
  example
  PRIVATE straw_project_options
  PUBLIC straw)
