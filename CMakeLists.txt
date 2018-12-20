cmake_minimum_required(VERSION 3.11)
project(CA_road)

set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
find_package(GSL REQUIRED)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_PREFIX_PATH  "~/libraries/")

find_library(ZOLTAN
        NAMES zoltan libzoltan
        HINTS /usr/local ~/local/lib ${CMAKE_PREFIX_PATH}/Zoltan/lib
        PATHS ${ZOLTAN_LIBRARY_PATH})

find_path(ZOLTAN_INCLUDES
        NAMES zoltan zoltan.h
        HINTS ${CMAKE_PREFIX_PATH}/Zoltan/include)

add_executable(CA_road
        include/vehicle.hpp
        include/ca_cell.hpp
        include/driving_direction.hpp
        include/rules.hpp
        src/main.cpp include/communication.hpp)
message(${CMAKE_SOURCE_DIR})
target_include_directories(CA_road INTERFACE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(CA_road PRIVATE ${ZOLTAN_INCLUDES})
target_link_libraries(CA_road ${ZOLTAN})
