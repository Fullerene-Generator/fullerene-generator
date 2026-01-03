#include <validation.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/internal/catch_compiler_capabilities.hpp>
#include <catch2/internal/catch_preprocessor_internal_stringify.hpp>
#include <catch2/internal/catch_test_registry.hpp>
#include <catch2/internal/catch_assertion_handler.hpp>
#include <catch2/internal/catch_decomposer.hpp>
#include <catch2/internal/catch_result_type.hpp>
#include <catch2/internal/catch_test_macro_impl.hpp>

#include <expansions/l_expansion.h>
#include <expansions/l_reduction.h>
#include <fullerene/construct.h>

#include <vector>
#include <string>

static std::vector<std::vector<unsigned int>> snapshot_rotations(const dual_fullerene& G)
{
    const std::size_t n = G.total_nodes();
    std::vector<std::vector<unsigned int>> rot(n);

    for (std::size_t id = 0; id < n; ++id) {
        auto v = G.get_node(static_cast<unsigned int>(id));
        rot[id].reserve(v->degree());
        for (std::size_t k = 0; k < v->degree(); ++k) {
            rot[id].push_back(v->neighbor_at(k)->id());
        }
    }

    return rot;
}

static void require_same_graph(const dual_fullerene& A, const dual_fullerene& B)
{
    REQUIRE(A.total_nodes() == B.total_nodes());

    const auto ra = snapshot_rotations(A);
    const auto rb = snapshot_rotations(B);

    REQUIRE(ra.size() == rb.size());
    for (std::size_t id = 0; id < ra.size(); ++id) {
        REQUIRE(ra[id].size() == rb[id].size());
        REQUIRE(ra[id] == rb[id]);
    }
}

static l_reduction matching_reduction_from_candidate(const dual_fullerene& G, const l_expansion& expan)
{
    l_reduction r;;
    r.first_edge = expan.inverse_first_edge();
    r.use_next = expan.candidate().clockwise;
    r.size = expan.candidate().length + 1;
    r.second_edge = expan.inverse_second_edge();

    return r;
}

static void run_one_case(dual_fullerene(*make)(), int i, const std::string& tag)
{
    auto G = make();
    validate_dual_fullerene(G);

    auto before = G;
   

    auto exps = find_l_expansions(G, i);
    INFO(tag);
    INFO("i = " << i);
    REQUIRE_FALSE(exps.empty());

    auto* le = dynamic_cast<l_expansion*>(exps.front().get());
    REQUIRE(le != nullptr);
    REQUIRE(le->validate());

    const auto c = le->candidate();

    le->apply();
    validate_dual_fullerene(G);

   

    const auto red = matching_reduction_from_candidate(G, *le);
    red.apply(G, c);

    validate_dual_fullerene(G);

    require_same_graph(G, before);
}

TEST_CASE("L reductions undo L0 expansions exactly (C20)") {
    run_one_case(&create_c20_fullerene, 0, "C20");
}

TEST_CASE("L reductions undo L1 expansions exactly (C20)") {
    run_one_case(&create_c20_fullerene, 1, "C20");
}


TEST_CASE("L reductions undo L expansions exactly (C28)") {
    run_one_case(&create_c28_fullerene, 0, "C28");
    run_one_case(&create_c28_fullerene, 1, "C28");
    run_one_case(&create_c28_fullerene, 3, "C28");
}
