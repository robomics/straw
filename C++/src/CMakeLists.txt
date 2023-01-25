find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)

add_library(straw INTERFACE)

target_sources(
    straw
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/straw_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_stream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/matrix_zoom_data_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/straw/straw.h
              ${CMAKE_CURRENT_SOURCE_DIR}/include/straw/internal/common.h)

target_include_directories(straw INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(straw INTERFACE straw_project_options CURL::libcurl ZLIB::ZLIB)
