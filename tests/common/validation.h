#ifndef VALIDATION_H
#define VALIDATION_H

#include <map>
#include <catch2/catch_test_macros.hpp>
#include <catch2/internal/catch_assertion_handler.hpp>
#include <catch2/internal/catch_decomposer.hpp>
#include <catch2/internal/catch_test_macro_impl.hpp>
#include <catch2/internal/catch_preprocessor_internal_stringify.hpp>
#include <fullerene/construct.h>
#include <fullerene/dual_fullerene.h>

inline void validate_dual_fullerene(const dual_fullerene& f) {
    REQUIRE(f.get_nodes_5().size() == 12);

    for (auto& u: f.get_nodes_5()) {
        REQUIRE(u->degree() == 5);

        for (auto& weak_v : u->neighbors()) {
            auto v = weak_v.lock();
            REQUIRE(v);

            size_t count = 0;
            for (auto& weak_v2 : u->neighbors())
                if (weak_v2.lock() == v) count++;

            REQUIRE(count == 1);
        }
    }

    for (auto& u: f.get_nodes_6()) {
        REQUIRE(u->degree() == 6);

        for (auto& weak_v : u->neighbors()) {
            auto v = weak_v.lock();
            REQUIRE(v);

            size_t count = 0;
            for (auto& weak_v2 : u->neighbors())
                if (weak_v2.lock() == v) count++;

            REQUIRE(count == 1);
        }
    }
}

inline void validate_fullerene(const fullerene& f, const dual_fullerene& d) {
    const auto& faces = f.get_adjacency();
    const std::size_t V = d.total_nodes();
    const std::size_t E = (5 * 12 + 6 * (V - 12)) / 2;
    const std::size_t F = E - V + 2;

    REQUIRE(faces.size() == F);

    auto normalize = [](unsigned a, unsigned b) {
        return std::minmax(a, b);
    };

    std::map<std::pair<unsigned, unsigned>, unsigned> edge_count;

    for (unsigned u = 0; u < faces.size(); ++u) {
        const auto& adj = faces[u];

        for (unsigned v : adj) {
            REQUIRE(std::find(faces[v].begin(), faces[v].end(), u) != faces[v].end());
        }

        const unsigned a = adj[0];
        const unsigned b = adj[1];
        const unsigned c = adj[2];

        REQUIRE(a != u);
        REQUIRE(b != u);
        REQUIRE(c != u);

        REQUIRE(a < F);
        REQUIRE(b < F);
        REQUIRE(c < F);

        REQUIRE(a != b);
        REQUIRE(b != c);
        REQUIRE(c != a);

        auto e1 = normalize(u, a);
        auto e2 = normalize(u, b);
        auto e3 = normalize(u, c);

        edge_count[e1]++;
        edge_count[e2]++;
        edge_count[e3]++;
    }
}

#endif //VALIDATION_H
