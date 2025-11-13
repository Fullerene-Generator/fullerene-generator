cmake_minimum_required(VERSION 3.30)
project(fullerene_generator)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES src/base_node.cpp
        src/directed_edge.cpp
        src/dual_fullerene.cpp
        src/main.cpp
)

add_executable(fullerene_generator ${SOURCES})

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)