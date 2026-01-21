#include <expansions/base_reduction.h>
#include <expansions/signature_state.h>
#include <expansions/l_reduction.h>
#include <expansions/b_reduction.h>

#include <cstdint>
#include <memory>
#include <vector>
#include <iostream>
#include <typeinfo>




static std::uint32_t edge_neighborhood_code(const directed_edge& e, bool use_next)
{
    auto e0 = e;
    std::uint32_t code = 0;

    auto pack = [&](int deg) {
        code <<= 1;
        if (deg == 6) code |= 1u;
        };

    for (int k = 0; k < 5; k++) {
        pack(static_cast<int>(e0.to()->degree()));
        e0 = use_next ? e0.next_around() : e0.prev_around();
    }
    return code;
}

static std::uint32_t path_neighborhood_code(const base_reduction& r, int edges_len = 7)
{
    std::uint32_t code = 0;
    if (edges_len <= 0) return code;

    auto e = r.first_edge;

    for (int step = 0; step < edges_len; ++step) {
        e = r.use_next ? e.right_turn(3) : e.left_turn(3);

        auto ep = e.prev_around();
        auto en = e.next_around();

        int dp = static_cast<int>(ep.to()->degree());
        int dn = static_cast<int>(en.to()->degree());

        std::uint32_t bits = 0;
        if (r.use_next) {
            if (dp == 6) bits |= 1u;
            if (dn == 6) bits |= 2u;
        }
        else {
            if (dn == 6) bits |= 1u;
            if (dp == 6) bits |= 2u;
        }

        code <<= 2;
        code |= bits;
    }

    return code;
}

std::uint32_t base_reduction::x2_code() const
{
    return edge_neighborhood_code(first_edge, use_next);
}

std::uint32_t base_reduction::x3_code() const
{
    return edge_neighborhood_code(second_edge, use_next);
}

std::uint32_t base_reduction::x4_code() const
{
    return path_neighborhood_code(*this, 7);
}

void base_reduction::fill_signature_candidate(expansion_candidate& out) const
{
    out.start = first_edge;
    out.clockwise = use_next;
    out.path.clear();
    out.parallel_path.clear();
}

bool base_reduction::is_canonical(const dual_fullerene& G, int min_x0, int l1, int l2) const
{
    const int ref_x0 = x0();

    if (ref_x0 > 2 && !G.is_ipr()) {
        return false;
    }

    for (int s = min_x0; s < ref_x0; ++s) {
        auto smaller = find_all_reductions(G, s, -1, -1, true, -1, -1);
        if (!smaller.empty()) {
            return false;
        }
    }
    

    auto candidates = find_all_reductions(G, ref_x0, first_edge.from->id(), first_edge.index, use_next, l1, l2);
    if (candidates.empty()) return true;
    
    const int ref_x1 = x1();
    {
        std::vector<std::unique_ptr<base_reduction>> next;
        next.reserve(candidates.size());
        for (auto& r : candidates) {
            int v = r->x1();
            if (v < ref_x1) {
                return false;
            }
            if (v == ref_x1) next.push_back(std::move(r));
        }
        candidates.swap(next);
    }

    if (candidates.empty()) return true;
    
    const std::uint32_t ref_x2 = x2_code();
    {
        std::vector<std::unique_ptr<base_reduction>> next;
        next.reserve(candidates.size());
        for (auto& r : candidates) {
            std::uint32_t v = r->x2_code();
            if (v < ref_x2) {
                return false;
            }
            if (v == ref_x2) next.push_back(std::move(r));
        }
        candidates.swap(next);
    }
    if (candidates.empty()) return true;
    

    const std::uint32_t ref_x3 = x3_code();
    {
        std::vector<std::unique_ptr<base_reduction>> next;
        next.reserve(candidates.size());
        for (auto& r : candidates) {
            std::uint32_t v = r->x3_code();
            if (v < ref_x3) return false;
            if (v == ref_x3) next.push_back(std::move(r));
        }
        candidates.swap(next);
    }
    if (candidates.empty()) return true;

    const std::uint32_t ref_x4 = x4_code();
    {
        std::vector<std::unique_ptr<base_reduction>> next;
        next.reserve(candidates.size());
        for (auto& r : candidates) {
            std::uint32_t v = r->x4_code();
            if (v < ref_x4) return false;
            if (v == ref_x4) next.push_back(std::move(r));
        }
        candidates.swap(next);
    }
    if (candidates.empty()) return true;
    
    expansion_candidate ref_cand;
    fill_signature_candidate(ref_cand);
    signature_state ref_state(G, ref_cand);

    std::vector<expansion_candidate> cand_data;
    cand_data.resize(candidates.size());
    std::vector<signature_state> states;
    states.reserve(candidates.size());

    for (std::size_t i = 0; i < candidates.size(); ++i) {
        candidates[i]->fill_signature_candidate(cand_data[i]);
        states.emplace_back(G, cand_data[i]);
    }

    std::vector<bool> alive(states.size(), true);
    std::size_t alive_count = states.size();
    std::size_t prefix_len = 0;

    while (alive_count > 0) {
        if (!ref_state.finished()) ref_state.extend_step();
        for (std::size_t i = 0; i < states.size(); ++i) {
            if (!alive[i]) continue;
            if (!states[i].finished()) states[i].extend_step();
        }

        const auto& ref_sig = ref_state.signature();
        std::size_t ref_len = ref_sig.size();

        bool any_progress = false;

        for (std::size_t i = 0; i < states.size(); ++i) {
            if (!alive[i]) continue;

            const auto& sig = states[i].signature();
            std::size_t cand_len = sig.size();

            std::size_t j = prefix_len;
            if (j >= ref_len && j >= cand_len) continue;

            any_progress = true;

            for (; j < ref_len; ++j) {
                int vr = ref_sig[j];
                int vc = sig[j];

                if (vc < vr) return false;
                if (vc > vr) {
                    alive[i] = false;
                    --alive_count;
                    break;
                }
            }
        }

        prefix_len = ref_len;
        if (!any_progress) break;
    }

    return true;
}

std::vector<std::unique_ptr<base_reduction>>
find_all_reductions(const dual_fullerene& G, int x0, int skip_pent, int skip_index, bool skip_clockwise, int skip_l1, int skip_l2)
{
    std::vector<std::unique_ptr<base_reduction>> out;

    const int l_param = x0;
    if (l_param >= 1) {
        std::vector<l_reduction> ls;
        if (skip_l1 >= 0) {
            ls = find_l_reductions(G, l_param, -1, -1, skip_clockwise);
        }
        else {
            ls = find_l_reductions(G, l_param, skip_pent, skip_index, skip_clockwise);
        }
        
        out.reserve(ls.size());
        for (auto& r : ls) {
            out.push_back(std::make_unique<l_reduction>(std::move(r)));
        }
    }
    const int b_sum = x0 - 2;
    for (int b = 0; b <= b_sum; b++) {
        std::vector<b_reduction> bs;
        if (b == skip_l1) {
            bs = find_b_reductions(G, b, b_sum - b, skip_pent, skip_index, skip_clockwise);
        }
        else {
            bs = find_b_reductions(G, b, b_sum - b, -1, -1, true);
        }
        
        out.reserve(out.size() + bs.size());
        for (auto& r : bs) {
            out.push_back(std::make_unique<b_reduction>(std::move(r)));
        }
    }

    return out;
}


int limit_by_reduction_distances(const dual_fullerene& G, int cur_best) {
    constexpr int N = 12;
    constexpr uint16_t FULL = (1u << N) - 1u; 
    std::array<uint8_t, 1u << N> seenPents{};
    std::array<uint8_t, 1u << N> seenPairUnion{};

    std::vector<uint16_t> pents;     pents.reserve(66);  
    std::vector<uint16_t> pairList;  pairList.reserve(495); // distinct 4-bit masks

    int distinctEdges = 0;

    auto reds = find_all_reductions(G, 2);

    for (const auto& r : reds) {
        int a = r->first_edge.from->id(), b = r->second_edge.from->id();
        uint16_t e = (uint16_t)((1u << a) | (1u << b));

        if (seenPents[e]) continue;
        seenPents[e] = 1;
        distinctEdges++;

        //check if it forms new disjoint triple

        for (uint16_t u : pairList) {
            if ((u & e) == 0) {
                return 1;
            }
        }

        // Update disjoint-pair unions
        for (uint16_t prev : pents) {
            if ((prev & e) == 0) {
                uint16_t u = (uint16_t)(prev | e); // 4-bit union
                if (!seenPairUnion[u]) {
                    seenPairUnion[u] = 1;
                    pairList.push_back(u);
                }
            }
        }

        pents.push_back(e);
    }

    if (distinctEdges >= 2) return 2;
    return cur_best;
}
