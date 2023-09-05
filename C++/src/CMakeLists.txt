find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)

option(STRAW_BUILD_EXECUTABLE "Build straw executable" ON)

# Add straw the library
add_library(straw)
add_library(straw::straw ALIAS straw)

target_sources(straw PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/straw.cpp")
target_include_directories(straw PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:include>")

target_link_libraries(straw PUBLIC CURL::libcurl ZLIB::ZLIB m)


if (STRAW_BUILD_EXECUTABLE)
    # Add straw the executable
    add_executable(straw_main)
    set_target_properties(straw_main PROPERTIES OUTPUT_NAME "straw")

    target_sources(straw_main PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
    target_link_libraries(straw_main PUBLIC straw::straw)
endif()
