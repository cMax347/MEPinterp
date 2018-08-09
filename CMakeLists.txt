#configure for fotran use
cmake_minimum_required (VERSION 2.6)
enable_language (Fortran)

#project name
project(MEPinterp)
#version info
set (cmakeTest_VERSION_MAJOR 1)
set (cmakeTest_VERSION_MINOR 0)


#set some directions
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/archive)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

#use all files provided in src, i.e. my .f90 files
file(GLOB_RECURSE sources src/*.f90 src/*.F90 )



message("detected system ${CMAKE_SYSTEM_NAME}")



#for LINUX use intel compiler with performance options
if (UNIX AND NOT APPLE)
	message("will use LINUX (intel) compiler options")
  set(CMAKE_Fortran_FLAGS "-mkl -O3 -qopenmp")
  set(COVERAGE_LINK_FLAGS "-mkl -O3 -qopenmp")
    #do something
endif (UNIX AND NOT APPLE)

#for APPLE use a debugging flag scheme
if (APPLE)
  message("will use MAC compiler options")
  set(CMAKE_Fortran_FLAGS "-Og -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -openmp")
  set(COVERAGE_LINK_FLAGS "-Og -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -openmp")
endif (APPLE)


find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  add_executable(main.exe ${sources})
  target_link_libraries(main.exe ${LAPACK_LIBRARIES})
  message("Lapack link library is ${LAPACK_LIBRARIES}")
else(LAPACK_FOUND)
    message("no lapack library found (mandatory)")
endif(LAPACK_FOUND)

find_package(MPI )
if(MPI_FOUND)
  target_link_libraries(main.exe ${MPI_LIBRARIES})
  message("mpi link library is ${MPI_LIBRARIES}")
else(MPI_FOUND)
  message("no mpi library found")
endif(MPI_FOUND)



message ("commpiler prepared")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")




#add thirdparty libraries
include_directories("${PROJECT_SOURCE_DIR}/thirdparty")
add_subdirectory(thirdparty)
target_link_libraries(main.exe config_fortran ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/wannier90/libwannier.a)
