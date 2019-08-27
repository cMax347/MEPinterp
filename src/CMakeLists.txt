
# add library with all src files, excluding main.f90
file(GLOB srcLst ./*.f90)
list(REMOVE_ITEM srcLst main.f90)
add_library(srcLib ${srcLst})
