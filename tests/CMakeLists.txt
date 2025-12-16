set(TEST_SOURCES
        test_fullerene.cpp
        test_expansions.cpp
        test_embeddings.cpp
)

add_executable(fullerene_tests ${TEST_SOURCES} "test_reductions.cpp")

target_link_libraries(fullerene_tests
        PRIVATE
        fullerene_core
        fullerene_expansions
        fullerene_generators
        embedding_lib
        Catch2::Catch2WithMain
)

target_include_directories(fullerene_tests
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        common
)

include(CTest)
include(Catch)
catch_discover_tests(fullerene_tests)
