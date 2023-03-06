cmake_minimum_required(VERSION 3.0.0)
project(SLAE VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

#enable_testing()

add_compile_options(-Wall -Werror)

add_subdirectory(source/csr_matrix)
add_subdirectory(source/tridiagonal_matrix)
add_subdirectory(source/iteration_methods)
add_subdirectory(TESTS)
