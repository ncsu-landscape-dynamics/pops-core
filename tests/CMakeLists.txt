# adds a .cpp file as a test
# takes one parameter which is a filename without an extension
function(add_pops_test NAME)
    # a test is an executable
    add_executable("${NAME}" "${NAME}.cpp")

    # make the PoPS library a dependency
    target_link_libraries(${NAME} pops)

    # Enable compiler warnings
    target_compile_options(${NAME} PRIVATE
         $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
              -Wall -Wextra -pedantic -Wfloat-conversion -Werror>
         $<$<CXX_COMPILER_ID:MSVC>:
              /W4>)
     if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
         add_compile_options(-Wdocumentation -Wdocumentation-pedantic -Wno-documentation-unknown-command)
     endif()

    # register a test
    add_test(NAME "${NAME}" COMMAND ${NAME})
endfunction()

add_pops_test(test_competency_table)
add_pops_test(test_date)
add_pops_test(test_deterministic)
add_pops_test(test_distributions)
add_pops_test(test_environment)
add_pops_test(test_generator_provider)
add_pops_test(test_model)
add_pops_test(test_mortality)
add_pops_test(test_movements)
add_pops_test(test_multi_host_model)
add_pops_test(test_network)
add_pops_test(test_network_helpers)
add_pops_test(test_network_kernel)
add_pops_test(test_overpopulation_movements)
add_pops_test(test_quarantine)
add_pops_test(test_random)
add_pops_test(test_raster)
add_pops_test(test_scheduling)
add_pops_test(test_simulation)
add_pops_test(test_simulation_kernels)
add_pops_test(test_soils)
add_pops_test(test_spread_rate)
add_pops_test(test_statistics)
add_pops_test(test_survival_rate)
add_pops_test(test_treatments)
