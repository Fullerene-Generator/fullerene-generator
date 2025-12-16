add_library(fullerene_generators
        f_expansion_generator.cpp
 "l_expansion_generator.cpp")

target_link_libraries(fullerene_generators PUBLIC fullerene_expansions)
target_include_directories(fullerene_generators PUBLIC ${PROJECT_SOURCE_DIR}/include)
