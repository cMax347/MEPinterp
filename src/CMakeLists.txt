
# add library with all src files, excluding main.f90
file(GLOB srcLst ./*.f90 ./*.F)
list(REMOVE_ITEM srcLst main.f90)
add_library(srcLib ${srcLst})


#link thirdparty libraries
target_link_libraries(srcLib m_config)

#build executable
add_executable(mepInterp main.f90)
target_link_libraries(mepInterp srcLib)
