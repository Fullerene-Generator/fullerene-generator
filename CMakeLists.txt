cmake_minimum_required(VERSION 3.30)
project(fullerene_generator CXX)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(vendor/eigen)

set(FULLERENE_SOURCES
        src/base_node.cpp
        src/directed_edge.cpp
        src/dual_fullerene.cpp
        src/fullerene.cpp
        src/l_expansion.cpp
        src/f_expansion.cpp
        src/f_expansion_generator.cpp
)

add_library(fullerene_lib STATIC ${FULLERENE_SOURCES})

target_include_directories(fullerene_lib
        PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fullerene_lib
        PUBLIC
        Eigen3::Eigen
)

add_executable(fullerene_generator src/main.cpp)

target_link_libraries(fullerene_generator
        PRIVATE
        fullerene_lib
)
