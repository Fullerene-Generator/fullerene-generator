cmake_minimum_required(VERSION 3.30)
project(fullerene_generator)

set(CMAKE_CXX_STANDARD 20)

set(SOURCES
    src/construct.cpp
    src/fullerene.cpp
    src/main.cpp
)

add_executable(fullerene_generator ${SOURCES} "include/expansions.h" "src/expansions.cpp")

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)