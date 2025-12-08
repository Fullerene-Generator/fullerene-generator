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
            if (!v)
                return false;
            bool found = false;
            for (auto& weak_u2 : v->neighbors()) {
                if (auto u2 = weak_u2.lock()) {
                    if (u2.get() == u.get()) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found)
                return false;
        }
    }

    for (auto& u: f.get_nodes_6()) {
        if (u->degree() != 6)
            return false;

        for (auto& weak_v : u->neighbors()) {
            auto v = weak_v.lock();
            if (!v)
                return false;
            bool found = false;
            for (auto& weak_u2 : v->neighbors()) {
                if (auto u2 = weak_u2.lock()) {
                    if (u2.get() == u.get()) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found)
                return false;
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

TEST_CASE()
