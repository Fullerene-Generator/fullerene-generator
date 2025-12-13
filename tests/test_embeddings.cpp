#include <validation.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/internal/catch_compiler_capabilities.hpp>
#include <catch2/internal/catch_preprocessor_internal_stringify.hpp>
#include <catch2/internal/catch_test_registry.hpp>
#include <catch2/internal/catch_assertion_handler.hpp>
#include <catch2/internal/catch_decomposer.hpp>
#include <embeddings/embedder.h>
#include <expansions/f_expansion.h>
#include <fullerene/construct.h>

constexpr double EPS = 1e-6;

// test embedder
TEST_CASE("Tutte embedding test for base fullerenes") {
    const auto d1 = create_c20_fullerene();
    const auto d2 = create_c28_fullerene();
    const auto d3 = create_c30_fullerene();

    const auto f1 = d1.to_primal();
    const auto f2 = d2.to_primal();
    const auto f3 = d3.to_primal();

    const auto g1 = graph {
        f1.get_adjacency(),
        f1.get_outer_face_nodes()
    };
    const auto g2 = graph {
        f2.get_adjacency(),
        f2.get_outer_face_nodes()
    };
    const auto g3 = graph {
        f3.get_adjacency(),
        f3.get_outer_face_nodes()
    };

    auto e1 = embedder::compute_tutte(g1);
    auto e2 = embedder::compute_tutte(g2);
    auto e3 = embedder::compute_tutte(g3);

    for (auto& coords: e1) {
        REQUIRE(hypot(coords[0], coords[1]) <= 1 + EPS);
    }
    for (auto& coords: e2) {
        REQUIRE(hypot(coords[0], coords[1]) <= 1 + EPS);
    }
    for (auto& coords: e3) {
        REQUIRE(hypot(coords[0], coords[1]) <= 1 + EPS);
    }
}
