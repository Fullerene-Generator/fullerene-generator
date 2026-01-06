add_library(embedding_lib
        embedder.cpp
        forces.cpp
)

target_include_directories(embedding_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(embedding_lib PUBLIC Eigen3::Eigen)
