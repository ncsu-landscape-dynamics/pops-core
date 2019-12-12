# adds a .cpp file as a test
# takes one parameter which is a filename without an extension
function(add_pops_test NAME)
    # a test is an executable
    add_executable("${NAME}" "${NAME}")

    # make the PoPS library a dependency
    target_link_libraries(${NAME} pops)

    # register a test
    add_test(NAME "${NAME}" COMMAND ${NAME})
endfunction()

add_pops_test(test_date)
#add_pops_test(test_mortality)
add_pops_test(test_raster)
add_pops_test(test_simulation)
add_pops_test(test_spread_rate)
add_pops_test(test_statistics)
add_pops_test(test_treatments)
