cmake_minimum_required(VERSION 3.20.0)

project(polar_project)

find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(polar polar.cpp)
target_link_libraries(polar Eigen3::Eigen)
target_compile_options(polar PRIVATE -ggdb)

add_executable(ballistics ballistics.cpp conversions.h caams.cpp caams.hpp)
target_link_libraries(ballistics Eigen3::Eigen)
target_compile_options(ballistics PRIVATE -O3 -ggdb)


