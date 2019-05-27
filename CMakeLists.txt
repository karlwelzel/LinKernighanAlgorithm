cmake_minimum_required(VERSION 3.13)
project(LinKerninghanAlgorithm)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -pedantic-errors")

add_executable(LinKerninghanAlgorithm main.cpp Tour.cpp Tour.h TsplibUtils.cpp TsplibUtils.h
        SimpleHeuristic.cpp SimpleHeuristic.h LinKernighanHeuristic.cpp LinKernighanHeuristic.h
        SignedPermutation.cpp SignedPermutation.h)