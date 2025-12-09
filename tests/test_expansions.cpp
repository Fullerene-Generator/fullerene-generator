#include <validation.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/internal/catch_compiler_capabilities.hpp>
#include <catch2/internal/catch_preprocessor_internal_stringify.hpp>
#include <catch2/internal/catch_test_registry.hpp>
#include <catch2/internal/catch_assertion_handler.hpp>
#include <catch2/internal/catch_decomposer.hpp>
#include <catch2/internal/catch_result_type.hpp>
#include <catch2/internal/catch_test_macro_impl.hpp>
#include <expansions/f_expansion.h>
#include <fullerene/construct.h>

constexpr int EXPANSION_LIMIT = 10;

// test f_expansion
TEST_CASE("Test dual fullerene correctness after performing f_expansion") {
    auto d = create_c30_fullerene();
    validate_dual_fullerene(d);

    for (auto i = 0; i < EXPANSION_LIMIT; i++) {
        auto expansion = f_expansion(d, d.get_nodes_5()[0]);

        REQUIRE(expansion.validate());

        validate_dual_fullerene(d);
    }
}

TEST_CASE("Test (primal) fullerene correctness after performing f_expansion") {
    auto d = create_c30_fullerene();
    validate_dual_fullerene(d);

    for (auto i = 0; i < EXPANSION_LIMIT; i++) {
        auto expansion = f_expansion(d, d.get_nodes_5()[0]);

        REQUIRE(expansion.validate());

        validate_fullerene(d.to_primal(), d);
    }
}

TEST_CASE("Test f_expansion validation") {
    auto d = create_c30_fullerene();
    const auto e1 = f_expansion(d, d.get_nodes_5()[0]);
    const auto e2 = f_expansion(d, d.get_nodes_5()[1]);

    REQUIRE(e1.validate());
    REQUIRE(!e2.validate());
}
