# cmake_minimum_required(VERSION <specify CMake version here>)
project(rungekutta)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES
    main.cpp
    Parameters.cpp
    Parameters.h
    RK4.cpp
    RK4.h)

add_executable(rungekutta ${SOURCE_FILES})