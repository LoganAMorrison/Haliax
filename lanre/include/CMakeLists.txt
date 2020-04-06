cmake_minimum_required(VERSION 3.11)
project(Lanre CXX)
add_subdirectory(lanre)

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

# Define two variables in order not to repeat ourselves.
set(LIBRARY_NAME lanre)

# Define the library
add_library(${LIBRARY_NAME} INTERFACE)

target_include_directories(${LIBRARY_NAME} INTERFACE /include)

target_link_libraries(${LIBRARY_NAME} PUBLIC)
