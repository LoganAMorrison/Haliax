set(CMAKE_CXX_STANDARD 11)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif ()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif ()

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(test_diffeq_stiff test_diffeq_stiff.cpp)
target_link_libraries(test_diffeq_stiff PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_diffeq_stiff PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_dual test_dual.cpp)
target_link_libraries(test_dual PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_dual PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_higgs_portal test_higgs_portal.cpp)
target_link_libraries(test_higgs_portal PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_higgs_portal PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_kinetic_mixing test_kinetic_mixing.cpp)
target_link_libraries(test_kinetic_mixing PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_kinetic_mixing PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_diffeq_non_stiff test_diffeq_non_stiff.cpp)
target_link_libraries(test_diffeq_non_stiff PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_diffeq_non_stiff PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_dm_models_constant_thermal_cross_section test_dm_models_constant_thermal_cross_section.cpp)
target_link_libraries(test_dm_models_constant_thermal_cross_section PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_dm_models_constant_thermal_cross_section PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_dm_models_darksun test_dm_models_darksun.cpp)
target_link_libraries(test_dm_models_darksun PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_dm_models_darksun PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_dm_models_dipole_dm test_dm_models_dipole_dm.cpp)
target_link_libraries(test_dm_models_dipole_dm PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_dm_models_dipole_dm PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_cosmology_thermodynamic_particle test_cosmology_thermodynamic_particle.cpp)
target_link_libraries(test_cosmology_thermodynamic_particle PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_cosmology_thermodynamic_particle PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_cosmology_standard_model test_cosmology_standard_model.cpp)
target_link_libraries(test_cosmology_standard_model PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_cosmology_standard_model PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_roots_bracket test_roots_bracket.cpp)
target_link_libraries(test_roots_bracket PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_roots_bracket PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_special_functions_lambertw test_special_functions_lambertw.cpp)
target_link_libraries(test_special_functions_lambertw PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_special_functions_lambertw PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_special_functions_besselk test_special_functions_besselk.cpp)
target_link_libraries(test_special_functions_besselk PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_special_functions_besselk PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_dm_models_compute_xi test_dm_models_compute_xi.cpp)
target_link_libraries(test_dm_models_compute_xi PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_dm_models_compute_xi PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_phase_space_rambo test_phase_space_rambo.cpp)
target_link_libraries(test_phase_space_rambo PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_phase_space_rambo PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_integrate_vegas test_integrate_vegas.cpp)
target_link_libraries(test_integrate_vegas PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_integrate_vegas PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_integrate_qags test_integrate_qags.cpp)
target_link_libraries(test_integrate_qags PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_integrate_qags PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_interpolate_univariate_spline test_interpolate_univariate_spline.cpp)
target_link_libraries(test_interpolate_univariate_spline PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_interpolate_univariate_spline PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)


add_executable(test_gamma_ray_spectra_neutral_pion test_gamma_ray_spectra_neutral_pion.cpp)
target_link_libraries(test_gamma_ray_spectra_neutral_pion PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_gamma_ray_spectra_neutral_pion PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_optimize_hybrd test_optimize_hybrd.cpp)
target_link_libraries(test_optimize_hybrd PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_optimize_hybrd PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

add_executable(test_optimize_hybrj test_optimize_hybrj.cpp)
target_link_libraries(test_optimize_hybrj PUBLIC
        lanre
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES})

set_target_properties(test_optimize_hybrj PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)