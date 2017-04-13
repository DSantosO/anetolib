cmake_minimum_required(VERSION 3.0.0)

project(aneto)



if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT CMAKE_VERBOSE_MAKEFILE)
	set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

include_directories(.)
include_directories(/usr/include/)
include_directories(/usr/include/eigen3/)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


FIND_PACKAGE(Boost COMPONENTS program_options system filesystem REQUIRED)
FIND_PACKAGE(MPFR REQUIRED )
#FIND_PACKAGE(Eigen REQUIRED )


set (MAKE_CURRENT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
message(STATUS ${CMAKE_CURRENT_LIST_DIR})

if(NOT CMAKE_CXX_COMPILER)
	set(CMAKE_CXX_COMPILER "/usr/bin/g++")
endif()

set(CMAKE_CXX_FLAGS "-std=c++11  -lmpfr -lgmp -lboost_system -lboost_filesystem -lboost_program_options -fext-numeric-literals")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -ggdb -O0 -Wall -Wextra -Wpedantic")

option(QUAD "Include quadruple precision in the testing." ON)
if (QUAD)
	message(STATUS "Quadruple precision is enable. If you have any problem run cmake with '-DQUAD=OFF' option.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lquadmath -DQUAD=1")
else()
	message(STATUS "Quadruple precision is disable and it will not be tested.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQUAD=0")
endif()

option(OMP "Compile with OpenMP." ON)
if (OMP)
	message(STATUS "OpenMP is enable. If you have any problem run cmake with '-DOMP=OFF' option.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DOMP=1")
else()
	message(STATUS "OpenMP is disable and it will not be tested.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMP=0")
endif()


# Doxygen doc compilation - Use with 'make doc'
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file("${CMAKE_CURRENT_LIST_DIR}/Doxyfile" "${CMAKE_CURRENT_LIST_DIR}/doc/Doxyfile" @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_LIST_DIR}/doc/Doxyfile"
WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/doc/"
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)


# Build tests
file(GLOB TESTS_SRC  "tests/src/*.cpp")
add_executable(tests ${TESTS_SRC})

