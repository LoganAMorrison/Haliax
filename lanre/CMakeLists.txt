cmake_minimum_required(VERSION 3.11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../build)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/build)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../lib")

project("Lanre" VERSION 0.1.0
        DESCRIPTION "Collection of numerical and high-energy physics routines for c++"
        LANGUAGES CXX)

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_BINARY_DIR}/include
        ${Boost_INCLUDE_DIRS}
)

set(CMAKE_MACOSX_RPATH ON)

add_subdirectory(include)
add_subdirectory(test)
