#include <cassert>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include "fullerene/dual_fullerene.h"
#include "fullerene/construct.h"
#include "expansions/l_expansion.h"
#include "expansions/l_signature_state.h"
#include "expansions/l_reduction.h"

using Signature = std::vector<int>;

struct SignatureLess {
    bool operator()(const Signature& a, const Signature& b) const {
        if (a.size() != b.size()) {
            return a.size() < b.size();
        }
        return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    }
};

static bool is_valid_dual_fullerene(const dual_fullerene& G) {
    bool ok = true;

    assert(G.get_nodes_5().size() == 12);

    G.for_each_node([&](const std::shared_ptr<base_node>& node) {
        if (node->degree() != node->expected_degree()) {
            ok = false;
            return;
        }

        for (std::size_t k = 0; k < node->degree(); ++k) {
            auto neigh = node->neighbor_at(k);
            if (!neigh->is_neighbor_of(node)) {
                std::cout << "neigbor missing " << neigh->id() << "degree: " << neigh->degree()
                    << " is not neighbor of " << node->id() << " degree: " << node->degree() << '\n';
                ok = false;
                return;
            }
        }
        });

    return ok;
}

static std::size_t count_distinct_signatures(const dual_fullerene& G, int i) {
    auto candidates = find_L_candidates(const_cast<dual_fullerene&>(G), i);

    std::vector<LSignatureState> states;
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

    std::set<Signature, SignatureLess> uniq;
    for (const auto& st : states) {
        uniq.insert(st.signature());
    }
    return uniq.size();
}

static void test_c20_L0_counts() {
    dual_fullerene G = create_c20_fullerene();

    auto candidates = find_L_candidates(G, 0);
    std::cout << "C20: L0 candidates = " << candidates.size() << "\n";

    assert(candidates.size() == 120);

    dual_fullerene G2 = create_c20_fullerene();
    auto expansions = find_L_expansions(G2, 0);
    std::cout << "C20: unique L0 expansions = " << expansions.size() << "\n";

    assert(!expansions.empty());
    assert(expansions.size() <= candidates.size());
}

static void test_c20_L0_grouping_matches_signatures() {
    dual_fullerene G = create_c20_fullerene();
    std::size_t sig_classes = count_distinct_signatures(G, 0);
    std::cout << "C20: L0 signature classes = " << sig_classes << "\n";

    dual_fullerene G2 = create_c20_fullerene();
    auto expansions = find_L_expansions(G2, 0);

    std::cout << "C20: expansions returned = " << expansions.size() << "\n";
    assert(expansions.size() == sig_classes);
}

static void test_c28_L0() {
    dual_fullerene G = create_c28_fullerene();
    auto candidates = find_L_candidates(G, 0);
    std::cout << "C28: L0 candidates = " << candidates.size() << "\n";
    assert(!candidates.empty());

    std::size_t sig_classes = count_distinct_signatures(G, 0);
    std::cout << "C28: L0 signature classes = " << sig_classes << "\n";

    dual_fullerene G2 = create_c28_fullerene();
    auto expansions = find_L_expansions(G2, 0);
    std::cout << "C28: expansions returned = " << expansions.size() << "\n";

    assert(expansions.size() <= candidates.size());
    assert(expansions.size() == sig_classes);
}

static void test_c30_L0() {
    dual_fullerene G = create_c30_fullerene();
    auto candidates = find_L_candidates(G, 0);
    std::cout << "C30: L0 candidates = " << candidates.size() << "\n";
    assert(!candidates.empty());

    std::size_t sig_classes = count_distinct_signatures(G, 0);
    std::cout << "C30: L0 signature classes = " << sig_classes << "\n";

    dual_fullerene G2 = create_c30_fullerene();
    auto expansions = find_L_expansions(G2, 0);
    std::cout << "C30: L0 expansions returned = " << expansions.size() << "\n";

    assert(expansions.size() <= candidates.size());
    assert(expansions.size() == sig_classes);
}

static void test_c30_Li_apply_valid(int i) {
    std::cout << "\nC30: testing L" << i << " apply/validity\n";

    dual_fullerene G_base = create_c30_fullerene();
    auto candidates = find_L_candidates(G_base, i);
    std::cout << "C30: L" << i << " candidates = " << candidates.size() << "\n";
    assert(!candidates.empty());

    for (const auto& c : candidates) {
        dual_fullerene G = create_c30_fullerene();
        L_expansion exp(G, c);

        if (!exp.validate()) {
            continue;
        }

        exp.apply();
        assert(is_valid_dual_fullerene(G));
    }
}

static void test_c30_Li_grouping_matches_signatures(int i) {
    std::cout << "\nC30: testing L" << i << " grouping/signatures\n";

    dual_fullerene G = create_c30_fullerene();
    auto candidates = find_L_candidates(G, i);
    std::cout << "C30: L" << i << " candidates = " << candidates.size() << "\n";
    assert(!candidates.empty());

    std::size_t sig_classes = count_distinct_signatures(G, i);
    std::cout << "C30: L" << i << " signature classes = " << sig_classes << "\n";

    dual_fullerene G2 = create_c30_fullerene();
    auto expansions = find_L_expansions(G2, i);
    std::cout << "C30: L" << i << " expansions returned = " << expansions.size() << "\n";

    assert(expansions.size() <= candidates.size());
    assert(expansions.size() == sig_classes);
}

static void print_L0_representatives(const dual_fullerene& G) {
    auto candidates = find_L_candidates(const_cast<dual_fullerene&>(G), 0);
    auto expansions = find_L_expansions(const_cast<dual_fullerene&>(G), 0);

    std::cout << "\n=== Represented L0 expansions ===\n";

    for (auto& exp_ptr : expansions) {
        auto* e = static_cast<L_expansion*>(exp_ptr.get());
        const auto& c = e->candidate();

        std::cout << "Candidate: "
            << " start=" << c.start.from->id()
            << " -> " << c.start.to()->id()
            << " use_next=" << c.use_next
            << "\n";

        LSignatureState st(G, c);
        while (!st.finished()) st.extend_step();

        const auto& sig = st.signature();
        std::cout << "Signature: ";
        for (int x : sig) std::cout << x << " ";
        std::cout << "\n\n";
    }
}

static void test_c20_L0_reductions() {
    std::cout << "\nC20: L0 reductions after one L0 expansion\n";

    dual_fullerene G = create_c20_fullerene();
    auto expansions = find_L_expansions(G, 0);
    std::cout << "C20: L0 expansions = " << expansions.size() << "\n";
    assert(!expansions.empty());

    auto* e = static_cast<L_expansion*>(expansions[0].get());
    std::cout << "expansion with first edge: " << e->candidate().start.from->id() << " to: " << e->candidate().start.to()->id() << '\n';
    std::cout << "use next: " << e->candidate().use_next << '\n';
    e->apply();
    assert(is_valid_dual_fullerene(G));

    int red_size = 1;
    auto reductions = find_L_reductions(G, red_size);
    std::cout << "C20: L0 reductions of size " << red_size
        << " found = " << reductions.size() << "\n";

    for (std::size_t i = 0; i < reductions.size(); ++i) {
        const auto& r = reductions[i];
        std::cout << "  [" << i << "] "
            << "first_edge " << r.first_edge.from->id()
            << " -> " << r.first_edge.to()->id()
            << ", second_edge " << r.second_edge.from->id()
            << " -> " << r.second_edge.to()->id()
            << ", use_next=" << r.use_next
            << ", size=" << r.size
            << "\n";
    }
}

int main() {
    try {
        bool run_expansion_tests = false;
        bool run_signature_tests = false;
        bool run_reduction_tests = true;

        if (run_expansion_tests) {
            test_c20_L0_counts();
            test_c30_Li_apply_valid(1);
            test_c30_Li_apply_valid(3);
        }

        if (run_signature_tests) {
            test_c20_L0_grouping_matches_signatures();
            test_c28_L0();
            test_c30_L0();
            test_c30_Li_grouping_matches_signatures(1);
            test_c30_Li_grouping_matches_signatures(3);

            print_L0_representatives(create_c20_fullerene());
            print_L0_representatives(create_c28_fullerene());
            print_L0_representatives(create_c30_fullerene());
        }

        if (run_reduction_tests) {
            test_c20_L0_reductions();
        }
    }
    catch (const std::exception& ex) {
        std::cerr << "Exception: " << ex.what() << "\n";
        return 1;
    }

    std::cout << "Selected L tests finished.\n";
    return 0;
}
