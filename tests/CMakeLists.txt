
#create a library from all tests
add_library(Helpers helpers.f90)
target_link_libraries(Helpers srcLib)


#add TEST-EXECUTABLES
#
#	#1
add_executable(testIO test_file_io.f90)
target_link_libraries(testIO Helpers)
#
#	#2
add_executable(testCrossP test_crossp.f90)
target_link_libraries(testCrossP Helpers)
#
#	#3
add_executable(testMatMath test_matrix_math.f90)
target_link_libraries(testMatMath Helpers)



add_test(NAME I/O           COMMAND testIO        )
add_test(NAME parameters    COMMAND testCrossP    )
add_test(NAME matrix_math   COMMAND testMatMath   )