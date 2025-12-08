#include <catch2/catch_test_macros.hpp>
#include <catch2/internal/catch_preprocessor_internal_stringify.hpp>
#include <catch2/internal/catch_assertion_handler.hpp>
#include <catch2/internal/catch_decomposer.hpp>
#include <fullerene/construct.h>
#include <fullerene/dual_fullerene.h>

bool validate_dual_fullerene(const dual_fullerene& f) {
    if (f.get_nodes_5().size() != 12)
        return false;

    for (auto& u: f.get_nodes_5()) {
        if (u->degree() != 5)
            return false;

        for (auto& weak_v : u->neighbors()) {
            auto v = weak_v.lock();
            REQUIRE(v);
            size_t count = 0;
            for (auto& weak_v2 : u->neighbors())
                if (weak_v2.lock() == v) count++;
            if (count != 1) return false;
        }
    }

    for (auto& u: f.get_nodes_6()) {
        if (u->degree() != 6)
            return false;

        for (auto& weak_v : u->neighbors()) {
            auto v = weak_v.lock();
            REQUIRE(v);
            size_t count = 0;
            for (auto& weak_v2 : u->neighbors())
                if (weak_v2.lock() == v) count++;
            if (count != 1) return false;
        }
    }

    return true;
}

TEST_CASE("Base dual fullerenes are structurally valid", "[dual_fullerene]") {
    const auto d1 = create_c20_fullerene();
    const auto d2 = create_c28_fullerene();
    const auto d3 = create_c30_fullerene();

    REQUIRE(validate_dual_fullerene(d1));
    REQUIRE(validate_dual_fullerene(d2));
    REQUIRE(validate_dual_fullerene(d3));
}

TEST_CASE("base_node basic functionality") {
    const auto n5 = node_5::create(1);
    const auto n6 = node_6::create(2);

    REQUIRE(n5->expected_degree() == 5);
    REQUIRE(n6->expected_degree() == 6);

    REQUIRE(n5->degree() == 0);
    REQUIRE(n6->degree() == 0);
}
