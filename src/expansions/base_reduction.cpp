#include <expansions/base_reduction.h>
#include <expansions/signature_state.h>
#include <expansions/l_reduction.h>

#include <cstdint>
#include <memory>
#include <vector>
#include <iostream>

static std::uint32_t edge_neighborhood_code(const directed_edge& e, bool use_next)
{
    //std::cout << "computing edge neighbourhood code\n";
    auto e0 = e;
    std::uint32_t code = 0;
   // std::cout << "now pack\n";

    auto pack = [&](int deg) {
        code <<= 1;
        if (deg == 6) code |= 1u;
        };
    //std::cout << "now loop\n";
    for (int k = 0; k < 5; k++) {
      //  std::cout << "k: " << k << " to degree: " << e0.to()->degree() << '\n';
        pack(static_cast<int>(e0.to()->degree()));
        e0 = use_next ? e0.next_around() : e0.prev_around();
    }
   // std::cout << "ready\n";
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

bool base_reduction::is_canonical(const dual_fullerene& G, int min_x0) const
{
   // std::cout << "in is canonical\n";
    //std::cout << "reduction first edge from: " << first_edge.from->id() << " to: " << first_edge.to()->id() << '\n';
    const int ref_x0 = x0();
    if (ref_x0 <= 0) return true;

    for (int s = min_x0; s < ref_x0; ++s) {
        auto smaller = find_all_reductions(G, s);
        if (!smaller.empty()) return false;
    }

    auto candidates = find_all_reductions(G, ref_x0);
   // std::cout << "found " << candidates.size() << " reduction candidates of size " << ref_x0 << '\n';
    if (candidates.empty()) return true;

    const int ref_x1 = x1();
    {
        std::vector<std::unique_ptr<base_reduction>> next;
        next.reserve(candidates.size());
        for (auto& r : candidates) {
            int v = r->x1();
            if (v < ref_x1) return false;
            if (v == ref_x1) next.push_back(std::move(r));
        }
        candidates.swap(next);
    }
    //std::cout << "after x1 stage candidates size: " << candidates.size() << '\n';
    if (candidates.empty()) return true;
    //std::cout << "now will compute x2\n";
    const std::uint32_t ref_x2 = x2_code();
    //std::cout << "x2 code: " << ref_x2 << '\n';
    {
        std::vector<std::unique_ptr<base_reduction>> next;
        next.reserve(candidates.size());
        for (auto& r : candidates) {
            std::uint32_t v = r->x2_code();
            if (v < ref_x2) return false;
            if (v == ref_x2) next.push_back(std::move(r));
        }
        candidates.swap(next);
    }
    //std::cout << "after x2 stage candidates size: " << candidates.size() << '\n';
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
    //std::cout << "after x3 stage candidates size: " << candidates.size() << '\n';
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
    //std::cout << "after x4 stage candidates size: " << candidates.size() << '\n';
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
find_all_reductions(const dual_fullerene& G, int x0)
{
    std::vector<std::unique_ptr<base_reduction>> out;

    const int l_param = x0;
    if (l_param >= 1) {
        auto ls = find_l_reductions(G, l_param);
        out.reserve(ls.size());
        for (auto& r : ls) {
            out.push_back(std::make_unique<l_reduction>(std::move(r)));
        }
    }

    return out;
}
