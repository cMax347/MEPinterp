#!/bin/bash

if [[ $PWD == *"MEPinterp/scripts"* ]]; then
	echo "it appears you try to run this script from within the scripts root folder (which would wipe all scripts!)"
	echo "if you are really keen to delete stuff in this folder, I suggest you do it manually"
else 
	echo "attempt to wipe build directory $PWD. remove ..."
	sleep 3
	#
	# remove cmake history 
	rm -rfv CMakeCache.txt CTestTestfile.cmake CMakeFiles ./*.py  cmake_install.cmake
	sleep 1 
	#
	# remove old build files & folders
	rm -rfv Makefile Testing archive modules src tests thirdparty
	sleep 1
	#
	# delete executables
	rm -rfv  bin/mepInterp bin/testIO bin/testMatMath bin/testkSpace
	sleep 1
	#
	echo "... all done. Wiped build directory $PWD"
fi 

