cmake_minimum_required(VERSION 3.0)
set(CMAKE_VERBOSE_MAKEFILE ON)
project(Haliax)

# Boost and its components
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

# Lanre: Core c++ library
add_subdirectory("lanre")

# Use PyBind11 to link Lanre to Haliax
#set(PYBIND11_CPP_STANDARD -std=c++1)
set(PYTHON_EXECUTABLE "/Users/loganmorrison/.pyenv/shims/python")
add_subdirectory("pybind11")


include_directories(
        "${CMAKE_SOURCE_DIR}/lanre/include"
        ${Boost_INCLUDE_DIRS}
)

set(SOURCE_DIR haliax)

pybind11_add_module(cosmology ${SOURCES} "${SOURCE_DIR}/cosmology.cpp")
target_link_libraries(
        cosmology PRIVATE
        lanre
        Eigen3::Eigen
        ${Boost_LIBRARIES}
)

pybind11_add_module(dm_models ${SOURCES} "${SOURCE_DIR}/dm_models.cpp") #target_link_libraries(models PRIVATE
target_link_libraries(
        dm_models PRIVATE
        lanre
        Eigen3::Eigen
        ${Boost_LIBRARIES}
)

pybind11_add_module(diffeq ${SOURCES} "${SOURCE_DIR}/diffeq.cpp") #target_link_libraries(models PRIVATE
target_link_libraries(
        diffeq PRIVATE
        lanre
        Eigen3::Eigen
        ${Boost_LIBRARIES}
)

pybind11_add_module(constants ${SOURCES} "${SOURCE_DIR}/constants.cpp") #target_link_libraries(models PRIVATE
target_link_libraries(
        constants PRIVATE
        lanre
        Eigen3::Eigen
        ${Boost_LIBRARIES}
)