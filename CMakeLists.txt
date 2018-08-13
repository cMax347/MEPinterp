#configure for fotran use
cmake_minimum_required (VERSION 2.6)
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

#use all files provided in src, i.e. my .f90 files
file(GLOB_RECURSE sources src/*.f90 src/*.F90 )



add_executable(mepInterp ${sources})



include_directories("${PROJECT_SOURCE_DIR}/thirdparty")
add_subdirectory(thirdparty)
#target_link_libraries(main.exe config_fortran wannier)
target_link_libraries(mepInterp config_fortran)








if(CMAKE_Fortran_COMPILER_IS_Intel)
  #use intel mkl for LAPACK

  set(CMAKE_Fortran_FLAGS "-mkl -O3 -openmp")
  set(COVERAGE_LINK_FLAGS "-mkl -O3 -openmp")




else(CMAKE_Fortran_COMPILER_IS_Intel) 
  #
  find_package(LAPACK REQUIRED)
  target_link_libraries(mepInterp ${LAPACK_LIBRARIES})
  #
  find_package(MPI REQUIRED)
  target_link_libraries(mepInterp ${MPI_LIBRARIES})
endif(CMAKE_Fortran_COMPILER_IS_Intel)


 





if(CMAKE_Fortran_COMPILER_IS_GNU) 
  set(CMAKE_Fortran_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
  set(COVERAGE_LINK_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
endif(CMAKE_Fortran_COMPILER_IS_GNU)

#for LINUX use intel compiler with performance options
#if (UNIX AND NOT APPLE)
# find_package(MPI)
#  if(MPI_FOUND)
#    target_link_libraries(mepInterp ${MPI_LIBRARIES})
#    message("mpi link library is ${MPI_LIBRARIES}")
#  else(MPI_FOUND)
#    message("no mpi library found")
#  endif(MPI_FOUND)
#
#
#
#  message("will use LINUX (intel) compiler options")
#  set(CMAKE_Fortran_FLAGS "-mkl -O3 -openmp")
#  set(COVERAGE_LINK_FLAGS "-mkl -O3 -openmp")
#    #do something
#endif (UNIX AND NOT APPLE)

#for APPLE use a debugging flag scheme
#if (APPLE)
#  message("will use MAC compiler options")
#  set(CMAKE_Fortran_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
#  set(COVERAGE_LINK_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
#endif (APPLE)







set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")




#add thirdparty libraries




 