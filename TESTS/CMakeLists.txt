cmake_minimum_required(VERSION 3.0.0)
project(SLAE VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)


find_package(GTest REQUIRED)

add_executable(test samost2.cpp)
add_test(test COMMAND samost2.cpp)
target_link_libraries(test gtest slae_lib)