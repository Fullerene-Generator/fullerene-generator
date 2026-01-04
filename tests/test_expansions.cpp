#include <set>
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
#include <expansions/l_expansion.h>
#include <expansions/l_reduction.h>
#include <expansions/signature_state.h>

#include "expansions/b_expansion.h"


constexpr int EXPANSION_LIMIT = 10;

// test f_expansion
TEST_CASE("Test dual fullerene correctness after performing f_expansion", "[f_expansion]") {
    auto d = create_c30_fullerene();
    validate_dual_fullerene(d);

    for (auto i = 0; i < EXPANSION_LIMIT; i++) {
        auto expansion = f_expansion(d, d.get_nodes_5()[0]);

        REQUIRE(expansion.validate());

        expansion.apply();
        validate_dual_fullerene(d);
    }
}

TEST_CASE("Test (primal) fullerene correctness after performing f_expansion", "[f_expansion]") {
    auto d = create_c30_fullerene();
    validate_dual_fullerene(d);

    for (auto i = 0; i < EXPANSION_LIMIT; i++) {
        auto expansion = f_expansion(d, d.get_nodes_5()[0]);

        REQUIRE(expansion.validate());

        expansion.apply();
        validate_fullerene(d.to_primal(), d);
    }
}

TEST_CASE("Test f_expansion validation", "[f_expansion]") {
    auto d = create_c30_fullerene();
    const auto e1 = f_expansion(d, d.get_nodes_5()[0]);
    const auto e2 = f_expansion(d, d.get_nodes_5()[1]);

    REQUIRE(e1.validate());
    REQUIRE(!e2.validate());
}

TEST_CASE("f_expansion increases fullerene size by expected amount", "[f_expansion]") {
    auto d = create_c30_fullerene();
    const auto initial = d.to_primal().get_adjacency().size();

    f_expansion e(d, d.get_nodes_5()[0]);
    REQUIRE(e.validate());

    e.apply();

    REQUIRE(d.to_primal().get_adjacency().size() == initial + 10);
}

TEST_CASE("No dangling neighbors after multiple F expansions", "[f_expansion]") {
    auto d = create_c30_fullerene();

    for (int i = 0; i < 5; ++i) {
        auto exp = f_expansion(d, d.get_nodes_5()[0]);
        REQUIRE(exp.validate());
        exp.apply();
    }

    d.for_each_node([&](const std::shared_ptr<base_node>& u) {
        for (auto& w : u->neighbors()) {
            REQUIRE(w.lock());
        }
    });
}


// test l_expansion
using signature = std::vector<int>;

struct SignatureLess {
    bool operator()(const signature& a, const signature& b) const {
        if (a.size() != b.size()) {
            return a.size() < b.size();
        }
        return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    }
};

static std::size_t count_distinct_signatures(const dual_fullerene& G, int i) {
    auto candidates = find_l_candidates(const_cast<dual_fullerene&>(G), i);

    std::vector<signature_state> states;
    states.reserve(candidates.size());
    for (const auto& c : candidates) {
        states.emplace_back(G, c);
    }

    bool any_not_finished = true;
    while (any_not_finished) {
        any_not_finished = false;
        for (auto& st : states) {
            if (!st.finished()) {
                st.extend_step();
                any_not_finished = true;
            }
        }
    }

    std::set<signature, SignatureLess> uniq;
    for (const auto& st : states) {
        uniq.insert(st.signature());
    }
    return uniq.size();
}

TEST_CASE("C20 L(0) candidate count is correct", "[l_expansion]") {
    dual_fullerene G = create_c20_fullerene();

    const auto candidates = find_l_candidates(G, 0);
    INFO("L0 candidates = " << candidates.size());

    REQUIRE(candidates.size() == 120);
}

TEST_CASE("C20 L(0) grouping matches signature classes", "[l_expansion]") {
    dual_fullerene G = create_c20_fullerene();
    std::size_t sig_classes = count_distinct_signatures(G, 0);

    dual_fullerene G2 = create_c20_fullerene();
    auto expansions = find_l_expansions(G2, 0);

    INFO("Signature classes = " << sig_classes);
    INFO("Returned expansions = " << expansions.size());

    REQUIRE(expansions.size() == sig_classes);
}

TEST_CASE("C30 L(i) expansions preserve dual fullerene validity", "[l_expansion]") {
    constexpr int i = 1;

    dual_fullerene base = create_c30_fullerene();
    auto candidates = find_l_candidates(base, i);

    REQUIRE_FALSE(candidates.empty());

    for (const auto& cand : candidates) {
        dual_fullerene G = create_c30_fullerene();
        l_expansion exp(G, cand);

        if (!exp.validate()) {
            continue;
        }

        exp.apply();
        validate_dual_fullerene(G);
    }
}

TEST_CASE("L-expansion inverse reduction exists and is canonical", "[L_reduction]") {
    constexpr int i = 0;

    dual_fullerene G0 = create_c20_fullerene();
    auto expansions = find_l_expansions(G0, i);

    REQUIRE_FALSE(expansions.empty());

    for (std::size_t idx = 0; idx < expansions.size(); ++idx) {
        const auto* base_e = static_cast<l_expansion*>(expansions[idx].get());
        const auto& cand = base_e->candidate();

        dual_fullerene G = create_c20_fullerene();
        l_expansion exp(G, cand);

        if (!exp.validate())
            continue;

        exp.apply();
        validate_dual_fullerene(G);

        int red_size = i + 1;
        auto reductions = find_l_reductions(G, red_size);

        bool match_found = false;

        for (const auto& r : reductions) {
            if (r.size != red_size) continue;
            if (r.use_next != cand.clockwise) continue;
            if (r.first_edge.from->id() != cand.start.from->id()) continue;
            if (r.second_edge.from->id() != static_cast<unsigned>(cand.parallel_path[red_size + 1])) continue;

            match_found = true;
            REQUIRE(r.is_canonical(G, red_size));
            break;
        }

        REQUIRE(match_found);
    }
}

// test b_expansion
TEST_CASE("C20 B(0, 0) candidate count is correct", "[b_expansion]") {
    dual_fullerene G = create_c20_fullerene();

    const auto candidates = find_b_candidates(G, 0, 0);
    INFO("B(0, 0) candidates = " << candidates.size());

    REQUIRE(candidates.size() == 120);
}

TEST_CASE("C20 B(0, 0) grouping matches signature classes", "[l_expansion]") {
    dual_fullerene G = create_c20_fullerene();
    std::size_t sig_classes = count_distinct_signatures(G, 0);

    dual_fullerene G2 = create_c20_fullerene();
    auto expansions = find_b_expansions(G2, 0, 0);

    INFO("Signature classes = " << sig_classes);
    INFO("Returned expansions = " << expansions.size());

    REQUIRE(expansions.size() == sig_classes);
}

TEST_CASE("C30 B(i, i) expansions preserve dual fullerene validity", "[l_expansion]") {
    constexpr int i = 1;

    dual_fullerene base = create_c30_fullerene();
    auto candidates = find_b_candidates(base, i, i);

    REQUIRE_FALSE(candidates.empty());

    for (const auto& cand : candidates) {
        dual_fullerene G = create_c30_fullerene();
        b_expansion exp(G, cand);

        if (!exp.validate()) {
            continue;
        }

        exp.apply();
        validate_dual_fullerene(G);
    }
}
