add_library(m_config config_fortran/m_config.f90)
add_library(m_npyff NPY-for-Fortran/src/npy.F90 NPY-for-Fortran/src/endian_swap.f90)
#link_directories(wannier90)

