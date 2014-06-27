cmake_minimum_required(VERSION 2.8.0) 
project(cloudtest Fortran)
find_package(PythonInterp REQUIRED)
enable_testing()

add_subdirectory(tests)
