#configure for fotran use
cmake_minimum_required (VERSION 2.6 FATAL_ERROR)
enable_language (Fortran)

#project name
project(MEPinterp)
#version info
set (cmakeTest_VERSION_MAJOR 1)
set (cmakeTest_VERSION_MINOR 0)


#set some directions
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/archive)


message(" found fortran mpi compiler id ${MPI_Fortran_COMPILER}")

if ("${MPI_Fortran_COMPILER_ID}" STREQUAL "GNU")
  message("found GNU MPI compiler")
endif()
#for LINUX use intel compiler with performance options
if (UNIX AND NOT APPLE)
  message("detected unix will assume intel compiler")
  set(CMAKE_Fortran_FLAGS "-mkl -fpp -O3 -openmp")
  set(COVERAGE_LINK_FLAGS "-mkl -fpp -O3 -openmp")
  add_definitions(-DUSE_MPI)
    #do something
endif (UNIX AND NOT APPLE)

#for APPLE use a debugging flag scheme
if (APPLE)
  message("will use MAC compiler options (requires Lapack)")
  find_package(LAPACK REQUIRED)
  find_package(MPI)
  if (MPI_FOUND)
    add_definitions(-DUSE_MPI)
    message("found MPI lib, will build parallel exe")
  endif(MPI_FOUND)

  set(CMAKE_Fortran_FLAGS "-Og -O0 -llapack -cpp -lblas -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
  set(COVERAGE_LINK_FLAGS "-Og -O0 -llapack -cpp -lblas -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
endif (APPLE)


if (WIN32)
  message("detected Windows operating system. You are on your own")
  find_package(LAPACK REQUIRED)
  find_package(MPI)
  if (MPI_FOUND)
    add_definitions(-DUSE_MPI)
    message("found MPI lib, will build parallel exe")
  endif(MPI_FOUND)

endif (WIN32)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")

#	ToDo's:
#
#	1) debug and release build
#	2) parallel and serial build
#	2) INTEL compiler matching
#


#GET THE thirdparty libraries
#include_directories(${PROJECT_SOURCE_DIR}/thirdparty)
#include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(thirdparty)
add_subdirectory(src)

#build tests
enable_testing()
add_subdirectory(tests)




