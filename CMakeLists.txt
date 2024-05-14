cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_COMPILER "g++-12")
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic -Winline -Wmain")

add_executable(exercise02 exercise02.cpp function.cpp output.cpp solver.cpp)
