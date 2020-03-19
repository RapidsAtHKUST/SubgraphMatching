cmake_minimum_required(VERSION 2.8)
project(SubgraphMatching)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -Wall -march=native -pthread")

add_subdirectory(graph)
add_subdirectory(matching)
add_subdirectory(utility)