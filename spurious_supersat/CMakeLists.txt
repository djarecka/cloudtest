cmake_minimum_required(VERSION 3.0)
project(pympdata CXX)

set(target "mpdata")

add_library(${target} SHARED ${target}.cpp)
set_target_properties(${target} PROPERTIES SUFFIX ".so")

# Python
find_package(PythonLibs)                                                        
if (NOT PYTHON_LIBRARIES)                                                       
  message(FATAL_ERROR "Python libraries not found")                                                                            
endif()                                                                         
target_include_directories(${target} PUBLIC ${PYTHON_INCLUDE_DIRS})             
target_link_libraries(${target} ${PYTHON_LIBRARIES}) 

# Boost Python & Timer
find_package(Boost COMPONENTS python timer system QUIET)                                     
if (NOT Boost_FOUND)                                                            
  message(FATAL_ERROR "Boost.Python not found")
endif()                                                                         
target_link_libraries(${target} ${Boost_LIBRARIES})  

# Blitz++
include(CheckCXXSourceCompiles)                                                 
check_cxx_source_compiles("                                                     
    #include <blitz/array.h>                                                    
    int main() {}                                                               
  " BLITZ_FOUND)                                                                
if (NOT BLITZ_FOUND)                                                            
  message(FATAL_ERROR "Blitz++ not found")
else()                                                                          
  target_link_libraries(${target} blitz)
endif() 

# libmpdata++
find_package(libmpdata++)                                                       
                                                                                
if (NOT libmpdataxx_FOUND)                                                      
  message(FATAL_ERROR "libmpdata++ not found")                                  
endif()                                                                         
                                                                                
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")                                       
  # make the Release flags the default                                          
  set(CMAKE_CXX_FLAGS ${libmpdataxx_CXX_FLAGS_RELEASE})                         
  # ignore CMake defaults                                                       
  set(CMAKE_CXX_FLAGS_RELEASE "")                                               
else()                                                                          
  set(CMAKE_CXX_FLAGS_DEBUG ${libmpdataxx_CXX_FLAGS_DEBUG})                     
endif() 
