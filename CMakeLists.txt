project(climateModeling)
cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_COMPILER "g++-12")
set(CMAKE_CXX_FLAGS "-std=c++20 -Wall -pedantic -Winline -Wmain")

add_executable(exercise03 exercise03.cpp function.cpp output.cpp solver.cpp)
