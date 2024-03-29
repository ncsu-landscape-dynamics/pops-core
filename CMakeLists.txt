# Works with X as minimum and tested through Y (X...Y)
cmake_minimum_required(VERSION 3.11...3.13)

project(pops-core
    VERSION 2.1
    DESCRIPTION "PoPS (Pest or Pathogen Spread) Model Core C++ library"
    LANGUAGES CXX)

# Only do these if this is the main project (and not as add_subdirectory)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # We need the following C++ standard
    set(CMAKE_CXX_STANDARD 17)

    # Disable compiler-specific extensions
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Useful for folders in some IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Enables testing for subdirectories and adds target test
    include(CTest)

    # Add documentation
    find_package(Doxygen)
    if(Doxygen_FOUND)
        set(DOXYGEN_EXTRACT_ALL YES)
        set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
        set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
        set(DOXYGEN_QT_AUTOBRIEF YES)
        set(DOXYGEN_SOURCE_BROWSER YES)
        set(DOXYGEN_TEMPLATE_RELATIONS YES)
        set(DOXYGEN_CALL_GRAPH YES)
        set(DOXYGEN_CALLER_GRAPH YES)
        set(DOXYGEN_INTERACTIVE_SVG YES)
        set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
        set(DOXYGEN_GENERATE_LATEX NO)

        # docs is the name of the target
        doxygen_add_docs(docs
            "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
            "${CMAKE_CURRENT_SOURCE_DIR}/CHANGELOG.md"
            "${CMAKE_CURRENT_SOURCE_DIR}/TECHNICALDEBT.md"
            "${CMAKE_CURRENT_SOURCE_DIR}/contributing_docs"
            "${PROJECT_SOURCE_DIR}/include"
        )
    else()
        message(STATUS "Doxygen not found, not building docs")
    endif()
endif()

add_library(pops INTERFACE)
target_include_directories(pops INTERFACE include/)
# Show files in IDEs
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    target_sources(pops INTERFACE
        include/pops/actions.hpp
        include/pops/environment.hpp
        include/pops/environment_interface.hpp
        include/pops/competency_table.hpp
        include/pops/config.hpp
        include/pops/host_pool.hpp
        include/pops/host_pool_interface.hpp
        include/pops/model_type.hpp
        include/pops/generator_provider.hpp
        include/pops/natural_anthropogenic_kernel.hpp
        include/pops/natural_kernel.hpp
        include/pops/anthropogenic_kernel.hpp
        include/pops/uniform_kernel.hpp
        include/pops/radial_kernel.hpp
        include/pops/treatments.hpp
        include/pops/raster.hpp
        include/pops/statistics.hpp
        include/pops/model.hpp
        include/pops/neighbor_kernel.hpp
        include/pops/deterministic_kernel.hpp
        include/pops/kernel.hpp
        include/pops/simulation.hpp
        include/pops/kernel_base.hpp
        include/pops/kernel_types.hpp
        include/pops/switch_kernel.hpp
        include/pops/spread_rate.hpp
        include/pops/date.hpp
        include/pops/scheduling.hpp
        include/pops/quarantine.hpp
        include/pops/pest_pool.hpp
        include/pops/power_law_kernel.hpp
        include/pops/hyperbolic_secant_kernel.hpp
        include/pops/logistic_kernel.hpp
        include/pops/exponential_power_kernel.hpp
        include/pops/exponential_kernel.hpp
        include/pops/cauchy_kernel.hpp
        include/pops/gamma_kernel.hpp
        include/pops/lognormal_kernel.hpp
        include/pops/network.hpp
        include/pops/network_kernel.hpp
        include/pops/normal_distribution_with_uniform_fallback.hpp
        include/pops/normal_kernel.hpp
        include/pops/multi_host_pool.hpp
        include/pops/pest_host_table.hpp
        include/pops/soils.hpp
        include/pops/weibull_kernel.hpp
    )
endif()

# Testing only available if this is the main app
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) AND BUILD_TESTING)
    add_definitions(-D POPS_TEST)  # TODO: remove the #ifdef from code
    add_subdirectory(tests)
endif()
