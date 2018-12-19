cmake_minimum_required(VERSION 3.11)
project(CA_road)

set(CMAKE_CXX_STANDARD 14)


find_package(MPI REQUIRED)
find_package(GSL REQUIRED)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_SOURCE_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
set(TEST_SOURCE_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
set(INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/includes)


set(CMAKE_PREFIX_PATH  "~/libraries/")

find_library(ZOLTAN
        NAMES zoltan libzoltan
        HINTS /usr/local ~/local/lib ${CMAKE_PREFIX_PATH}/Zoltan/lib
        PATHS ${ZOLTAN_LIBRARY_PATH})

find_path(ZOLTAN_INCLUDES
        NAMES zoltan zoltan.h
        HINTS ${CMAKE_PREFIX_PATH}/Zoltan/include
        )

add_executable(CA_road
        include/geometric_load_balancer.hpp
        include/unloading_model.hpp
        include/spatial_elements.hpp
        include/unloading_model.hpp
        include/utils.hpp
        include/zoltan_fn.hpp
        include/partitioner.hpp
        src/main.cpp)
target_link_libraries(CA_road ${ZOLTAN})
target_include_directories(CA_road PRIVATE ${ZOLTAN_INCLUDES})