cmake_minimum_required(VERSION 3.0.0)
project(SLAE VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)


find_package(GTest REQUIRED)

add_executable(test iteration_methods_test.cpp)
add_test(test COMMAND iteration_methods_test.cpp)
target_link_libraries(test gtest)