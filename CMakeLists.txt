cmake_minimum_required(VERSION 3.8)
project(LinKernighanAlgorithm)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -pedantic-errors")
#set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE "/bin/include-what-you-use;-Xiwyu;any")


add_executable(LinKernighanAlgorithm main.cpp
        Tour.cpp Tour.h
        TsplibUtils.cpp TsplibUtils.h
        LinKernighanHeuristic.cpp LinKernighanHeuristic.h
        SignedPermutation.cpp SignedPermutation.h
        PrimsAlgorithm.cpp PrimsAlgorithm.h
        AlphaDistances.cpp AlphaDistances.h)