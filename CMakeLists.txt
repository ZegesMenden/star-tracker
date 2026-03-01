cmake_minimum_required(VERSION 3.14)
project(star-tracker
    VERSION   0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Options ──────────────────────────────────────────────────────────────────

option(STARTRACKER_USE_FLOAT  "Use float instead of double"          ON)
option(STARTRACKER_BUILD_TESTS "Build unit tests"                    ON)

# ── Library ──────────────────────────────────────────────────────────────────

set(LIB_SOURCES
    src/attitude.cc
    src/camera.cc
    src/centroid.cc
    src/database.cc
    src/io.cc
    src/star.cc
    src/starid.cc
)

add_library(startracker ${LIB_SOURCES})

target_include_directories(startracker
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include/startracker>
)

target_compile_features(startracker PUBLIC cxx_std_17)

if(STARTRACKER_USE_FLOAT)
    target_compile_definitions(startracker PUBLIC STARTRACK_USE_FLOAT)
endif()

# ── Tools ────────────────────────────────────────────────────────────────────

option(STARTRACKER_BUILD_TOOLS "Build helper tools (generate-db, benchmark)" ON)

if(STARTRACKER_BUILD_TOOLS)
    add_executable(generate-db tools/generate-db.cc)
    target_link_libraries(generate-db PRIVATE startracker)

    add_executable(generate-stars tools/generate-stars.cc)
    target_include_directories(generate-stars PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/vendor
    )
    target_link_libraries(generate-stars PRIVATE startracker)

    add_executable(generate-catalog tools/generate-catalog.cc)
    target_link_libraries(generate-catalog PRIVATE startracker)

    add_executable(benchmark benchmark/benchmark.cc)
    target_include_directories(benchmark PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/vendor
    )
    target_link_libraries(benchmark PRIVATE startracker)
endif()

# ── Tests ────────────────────────────────────────────────────────────────────

if(STARTRACKER_BUILD_TESTS)
    enable_testing()

    # Catch2 single-header in vendor/
    add_library(Catch2 INTERFACE)
    target_include_directories(Catch2 INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/vendor
    )

    set(TEST_SOURCES
        test/test_main.cc
        test/test_attitude.cc
        test/test_camera.cc
        test/test_centroid.cc
        test/test_database.cc
        test/test_serialize.cc
        test/test_star.cc
        test/test_starid.cc
        test/test_vec.cc
    )

    add_executable(startracker_tests ${TEST_SOURCES})
    target_link_libraries(startracker_tests PRIVATE startracker Catch2)

    add_test(NAME startracker_tests COMMAND startracker_tests)
endif()

# ── Install ──────────────────────────────────────────────────────────────────

include(GNUInstallDirs)

install(TARGETS startracker
    EXPORT startrackerTargets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    FILES
        src/attitude.hpp
        src/camera.hpp
        src/centroid.hpp
        src/database.hpp
        src/io.hpp
        src/star.hpp
        src/starid.hpp
        src/stfloat.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/startracker
)

install(EXPORT startrackerTargets
    FILE        startrackerTargets.cmake
    NAMESPACE   startracker::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/startracker
)

# ── Package config ───────────────────────────────────────────────────────────

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/startrackerConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/startrackerConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/startracker
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/startrackerConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/startrackerConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/startrackerConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/startracker
)
