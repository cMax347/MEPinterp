#!/bin/bash

# remove cmake history 
rm -rf CMakeCache.txt CTestTestfile.cmake CMakeFiles ./*.py  cmake_install.cmake 

# remove old build files & folders
rm -rf Makefile Testing archive modules src tests thirdparty

# delete executables
rm -rf  bin/mepInterp bin/testIO bin/testMatMath bin/testkSpace