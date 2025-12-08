cmake_minimum_required(VERSION 3.30)
project(fullerene_generator)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(vendor/eigen)

set(SOURCES src/base_node.cpp
        src/directed_edge.cpp
        src/dual_fullerene.cpp
        src/main.cpp
        src/fullerene.cpp
        src/l_expansion.cpp
        src/f_expansion.cpp
        src/f_expansion_generator.cpp
        src/l_signature_state.cpp
)

add_executable(fullerene_generator ${SOURCES})

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fullerene_generator
        PRIVATE
        Eigen3::Eigen
)



add_executable(test_l0
    src/base_node.cpp
    src/directed_edge.cpp
    src/dual_fullerene.cpp
    src/fullerene.cpp
    src/l_expansion.cpp
    src/f_expansion.cpp
    src/f_expansion_generator.cpp
    src/l_signature_state.cpp
    tests/test_l0.cpp
)

target_include_directories(test_l0 PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_l0 PRIVATE
    Eigen3::Eigen
)