#configure for fotran use
cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)

#project name
project(MEPinterp)
#version info
set (cmakeTest_VERSION_MAJOR 1)
set (cmakeTest_VERSION_MINOR 0)


#set some directions
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)




#for LINUX use intel compiler with performance options
if (UNIX AND NOT APPLE)
  message("detected unix will assume intel compiler")
  set(CMAKE_Fortran_FLAGS "-mkl -O3 -openmp")
  set(COVERAGE_LINK_FLAGS "-mkl -O3 -openmp")
    #do something
endif (UNIX AND NOT APPLE)

#for APPLE use a debugging flag scheme
if (APPLE)
  message("will use MAC compiler options, search for Lapack and MPI library")
  find_package(LAPACK REQUIRED)
  find_package(MPI REQUIRED)

  set(CMAKE_Fortran_FLAGS "-Og -O0 -llapack -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
  set(COVERAGE_LINK_FLAGS "-Og -O0 -llapack -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
endif (APPLE)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")






#GET SOURCE FILES
file(GLOB test_src_lst  src/tests/*.f90)
file(GLOB mep_src_lst src/mepInterp/*.f90)
file(GLOB mep_main src/mepInterp/main.f90)
list(REMOVE_ITEM  mep_src_lst src/mepInterp/main.f90)

#GET THE thirdparty libraries
#include_directories(${PROJECT_SOURCE_DIR}/thirdparty)
add_subdirectory(thirdparty)



#EXECUTABLES
add_executable(mepInterp ${mep_main})
add_executable(testMEP ${test_src_lst})

#CREATE MEP LIBRARY
add_library(mep_src ${mep_src_lst})
target_link_libraries(mep_src config_fortran)
if(APPLE)
  target_link_libraries(mep_src ${LAPACK_LIBRARIES} ${MPI_LIBRARIES})
  target_link_libraries(mepInterp ${MPI_LIBRARIES})
  target_link_libraries(testMEP ${LAPACK_LIBRARIES})
endif(APPLE)

#LINK MEP-LIB TO EXECUTABLES
target_link_libraries(mepInterp mep_src)
target_link_libraries(testMEP mep_src)


