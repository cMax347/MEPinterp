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
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include_directories("${PROJECT_SOURCE_DIR}/thirdparty")


#for LINUX use intel compiler with performance options
if (UNIX AND NOT APPLE)
  message("detected unix will assume intel compiler")
  set(CMAKE_Fortran_FLAGS "-mkl -O3 -openmp")
  set(COVERAGE_LINK_FLAGS "-mkl -O3 -openmp")
    #do something
endif (UNIX AND NOT APPLE)

#for APPLE use a debugging flag scheme
if (APPLE)
  message("will use MAC compiler options")
  find_package(LAPACK REQUIRED)
  find_package(MPI REQUIRED)

  set(CMAKE_Fortran_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
  set(COVERAGE_LINK_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
endif (APPLE)






#use all files provided in src, i.e. my .f90 files
file(GLOB_RECURSE test_files       tests/*.f90 tests/*.F90 )

#put all sources(but main.f90) into a library
file(GLOB_RECURSE mep_src_lst src/*.f90)
list(REMOVE_ITEM  mep_src_lst src/main.f90)
add_library(mep_src ${mep_src_lst})



#     mepInterp   BUILD
add_executable(mepInterp src/main.f90)
add_subdirectory(thirdparty)
target_link_libraries(mepInterp mep_src config_fortran)
if(APPLE)
  target_link_libraries(mepInterp ${LAPACK_LIBRARIES} ${MPI_LIBRARIES})
endif(APPLE)


#     test      BUILD
add_executable(testMEP ${test_files})
target_link_libraries(testMEP mep_src)
if (APPLE)
  target_link_libraries(testMEP ${LAPACAK_LIBRARIES})
endif (APPLE)











#if(CMAKE_Fortran_COMPILER_IS_Intel)
#  #use intel mkl for LAPACK
#  set(CMAKE_Fortran_FLAGS "-mkl -O3 -openmp")
#  set(COVERAGE_LINK_FLAGS "-mkl -O3 -openmp")
#else(CMAKE_Fortran_COMPILER_IS_Intel) 
#  #
#  find_package(LAPACK REQUIRED)
#  target_link_libraries(mepInterp ${LAPACK_LIBRARIES})
#  #
#  find_package(MPI REQUIRED)
#  target_link_libraries(mepInterp ${MPI_LIBRARIES})
#endif(CMAKE_Fortran_COMPILER_IS_Intel)



#if(CMAKE_Fortran_COMPILER_IS_GNU) 
#  set(CMAKE_Fortran_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
#  set(COVERAGE_LINK_FLAGS "-Og -O0 -Wall -Wno-tabs -Wextra -Werror=line-truncation -pedantic -fimplicit-none -fcheck=all -fbacktrace ")
#endif(CMAKE_Fortran_COMPILER_IS_GNU)




set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")


#OLD:








 