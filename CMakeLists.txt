cmake_minimum_required(VERSION 3.12)
project(project)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(project main.cpp russians.h russians.cpp hirscherg.h hirschberg.cpp)