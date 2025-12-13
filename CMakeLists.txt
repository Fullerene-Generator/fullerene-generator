cmake_minimum_required(VERSION 3.30)
project(fullerene_generator CXX)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(vendor/eigen)
add_subdirectory(vendor/catch2)


add_subdirectory(src/fullerene)
add_subdirectory(src/expansions)
add_subdirectory(src/embeddings)
add_subdirectory(src/generators)
add_subdirectory(apps/generator)
add_subdirectory(apps/embedder)

enable_testing()
add_subdirectory(tests)

