
#create a library from all tests
file(GLOB testLst *.f90)
list(REMOVE_ITEM testLst main.f90)
list(REMOVE_ITEM testLst test_parameters.f90)
add_library(Tests ${testLst})


#link srcLib
target_link_libraries(Tests srcLib)

#create a test executable
add_executable(testExe main.f90)
target_link_libraries(testExe Tests)


add_executable(testCrossP test_crossp.f90)
target_link_libraries(testCrossP Tests)


