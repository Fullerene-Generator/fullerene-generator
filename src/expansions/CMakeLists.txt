add_library(fullerene_expansions
        l_expansion.cpp
        f_expansion.cpp
)

target_link_libraries(fullerene_expansions PUBLIC fullerene_core)
target_include_directories(fullerene_expansions PUBLIC ${PROJECT_SOURCE_DIR}/include)
