cmake_minimum_required(VERSION 3.1)
project(manifold)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 3rdparty dependenices
include_directories(3rdparty/Eigen)
include_directories(3rdparty/glm)
add_subdirectory(3rdparty/pybind11)

# Source code
include_directories(src/manifold)
add_subdirectory(src/manifold/manifold)
add_subdirectory(src/manifold/pybind)
