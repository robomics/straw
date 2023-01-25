find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)

add_library(straw)

target_sources(
    straw
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/straw.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_file.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_stream.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/matrix_zoom_data.cpp
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/straw/straw.h)

target_include_directories(straw PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
    straw
    PRIVATE straw_project_options
    PUBLIC CURL::libcurl ZLIB::ZLIB)
