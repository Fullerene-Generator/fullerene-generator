add_library(fullerene_expansions
        b_expansion.cpp
        l_expansion.cpp
        l_reduction.cpp
        base_reduction.cpp
        signature_state.cpp
        f_expansion.cpp
        base_expansion.cpp
)

target_link_libraries(fullerene_expansions PUBLIC fullerene_core)
target_include_directories(fullerene_expansions PUBLIC ${PROJECT_SOURCE_DIR}/include)
