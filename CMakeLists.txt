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

message("")
message("")
message("")


#
#
#-------------------------------------------------------------------------------------
#       GNU COMPILER (gfortran, mpifort)
#----------------------------------------------------------------------------------
#
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  message("will use GNU compiler flags")
  #
  #
  find_package(LAPACK REQUIRED)
  find_package(MPI)
  if (MPI_FOUND)
    add_definitions(-DUSE_MPI)    #preprocessor flag USE_MPI
    message("found MPI lib, set USE_MPI preprocessor macro")
  endif(MPI_FOUND)
  #
  set(CMAKE_Fortran_FLAGS "-Og -O0 -llapack -cpp -lblas -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
  set(COVERAGE_LINK_FLAGS "-Og -O0 -llapack -cpp -lblas -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
#
#-------------------------------------------------------------------------------------
#       INTEL COMPILER (mpiifort)
#----------------------------------------------------------------------------------
#
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  message("will use Intel compiler flags (assuming mpiifort)")
  add_definitions(-DUSE_MPI)
  message("INTEL compiler: set USE_MPI preprocessor macro")

  message("intel compiler version: ${CMAKE_Fortran_COMPILER_VERSION}")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "16.0.0.20161018")
    set(CMAKE_Fortran_FLAGS "-mkl -fpp -O3 -openmp")
    set(COVERAGE_LINK_FLAGS "-mkl -fpp -O3 -openmp")

    #set(CMAKE_Fortran_FLAGS "-mkl -fpp -O0 -g -traceback -check all -CB -check uninit -ftrapuv -debug all -openmp")
    #set(COVERAGE_LINK_FLAGS "-mkl -fpp -O0 -g -traceback -check all -CB -check uninit -ftrapuv -debug all -openmp")



    message("detected older intel compiler (openmp)")
  else()
    message("machine = ${MACHINE}")
    if("${MACHINE}" STREQUAL "jureca" )
      message("I think I am on jureca !!!1!")
    endif()
    set(CMAKE_Fortran_FLAGS "-mkl -fpp -O3 -qopenmp")
    set(COVERAGE_LINK_FLAGS "-mkl -fpp -O3 -qopenmp")
    message("detected newer intel compiler (qopenmp)")
  endif()
#
#-------------------------------------------------------------------------------------
#       other compilers
#----------------------------------------------------------------------------------
#
else ()
  message("unsupported compiler ID: ${CMAKE_Fortran_COMPILER_ID}, good luck!")
  find_package(LAPACK REQUIRED)
    if (MPI_FOUND)
    add_definitions(-DUSE_MPI)    #preprocessor flag USE_MPI
    message("found MPI lib, set USE_MPI preprocessor macro")
  endif(MPI_FOUND)
endif()






set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")

#	ToDo's:
#
#	1) debug and release build
#	2) parallel and serial build#


#GET THE thirdparty libraries
#include_directories(${PROJECT_SOURCE_DIR}/thirdparty)
#include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(thirdparty)
add_subdirectory(src)

#build tests
enable_testing()
add_subdirectory(tests)




