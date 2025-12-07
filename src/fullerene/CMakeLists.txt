add_library(fullerene_core
        base_node.cpp
        directed_edge.cpp
        dual_fullerene.cpp
        fullerene.cpp
)

target_include_directories(fullerene_core PUBLIC ${PROJECT_SOURCE_DIR}/include)
