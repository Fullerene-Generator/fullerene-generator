#include <expansions/l_reduction.h>
#include <expansions/l_signature_state.h>
#include <iostream>
#include <cstdint>
#include <vector>

std::vector<l_reduction>
find_L_reductions(const dual_fullerene& G,
    int size,
    int skip_pent_1,
    int skip_pent_2)
{
    std::vector<l_reduction> out;
    if (size < 1) {
        return out;
    }

    int path_len = size + 1;

    for (const auto& pent : G.get_nodes_5()) {
        auto start_node = pent;
        int deg = static_cast<int>(start_node->degree());

        for (int i = 0; i < deg; ++i) {
            directed_edge e0{ start_node, static_cast<std::size_t>(i) };

            for (bool use_next : { true, false }) {
                int outside_hex1 = use_next
                    ? static_cast<int>(e0.prev_around(2).to()->degree())
                    : static_cast<int>(e0.next_around(2).to()->degree());
                if (outside_hex1 != 6) {
                    continue;
                }

                auto e = e0;
                std::vector<int> path;
                path.reserve(path_len);

                bool ok = true;
                for (int step = 0; step < path_len; ++step) {
                    auto v = e.from;
                    unsigned int vid = v->id();

                    if (step == 0 || step == path_len - 1) {
                        if (v->type() != node_type::NODE_5) {
                            ok = false;
                            break;
                        }
                    }
                    else {
                        if (v->type() != node_type::NODE_6) {
                            ok = false;
                            break;
                        }
                    }

                    path.push_back(static_cast<int>(vid));

                    if (step < path_len - 1) {
                        if (use_next) {
                            e = e.right_turn(3);
                        }
                        else {
                            e = e.left_turn(3);
                        }
                    }
                }

                if (!ok) {
                    continue;
                }

                int outside_hex2 = use_next
                    ? static_cast<int>(e.next_around(1).to()->degree())
                    : static_cast<int>(e.prev_around(1).to()->degree());
                if (outside_hex2 != 6) {
                    continue;
                }

                int first_pent = path.front();
                int last_pent = path.back();

                if (skip_pent_1 >= 0 && skip_pent_2 >= 0) {
                    bool same_order = (first_pent == skip_pent_1 && last_pent == skip_pent_2);
                    bool swapped_order = (first_pent == skip_pent_2 && last_pent == skip_pent_1);
                    if (same_order || swapped_order) {
                        continue;
                    }
                }

                auto last_node = G.get_node(static_cast<unsigned int>(last_pent));
                auto prev_node = G.get_node(static_cast<unsigned int>(path[path.size() - 2]));
                directed_edge second_edge = last_node->get_edge(prev_node);

                l_reduction r;
                r.first_edge = e0;
                r.second_edge = second_edge;
                r.use_next = use_next;
                r.size = size;

                out.push_back(r);
            }
        }
    }

    return out;
}

static std::uint32_t edge_neighborhood_code(const directed_edge& e, bool use_next)
{
    auto e0 = e;
    std::uint32_t code = 0;
    auto pack = [&](int deg) {
        code <<= 1;
        if (deg == 6) {
            code |= 1u;
        }
        };

    for (int k = 0; k < 5; k++) {
        pack(e0.to()->degree());
        e0 = use_next ? e0.next_around() : e0.prev_around();
    }

    return code;
}


static std::uint32_t path_neighborhood_code(const l_reduction& r,
    int edges_len = 7)
{
    std::uint32_t code = 0;
    if (edges_len <= 0) {
        return code;
    }

    auto e = r.first_edge;

    for (int step = 0; step < edges_len; ++step) {
        if (r.use_next) {
            e = e.right_turn(3);
        }
        else {
            e = e.left_turn(3);
        }

        auto ep = e.prev_around();
        auto en = e.next_around();

        int dp = static_cast<int>(ep.to()->degree());
        int dn = static_cast<int>(en.to()->degree());

        std::uint32_t bits = 0;
        if (r.use_next) {
            if (dp == 6) {
                bits |= 1u;
            }
            if (dn == 6) {
                bits |= 2u;
            }
        }
        else {
            if (dn == 6) {
                bits |= 1u;
            }
            if (dp == 6) {
                bits |= 2u;
            }
        }
        

        code <<= 2;
        code |= bits;
    }

    return code;
}

static std::string to_binary(std::uint32_t x) {
    std::string s;
    for (int i = 31; i >= 0; --i) {
        s.push_back((x & (1u << i)) ? '1' : '0');
    }
    return s;
}

bool l_reduction::is_canonical(const dual_fullerene& G, int min_size) const
{
    int ref_size = size;
    if (ref_size <= 0) {
        return true;
    }

    for (int s = min_size; s < ref_size; ++s) {
        auto smaller = find_L_reductions(G, s);
        if (!smaller.empty()) {
            return false;
        }
    }

    int path_len = ref_size + 1;
    int ref_first_pent = static_cast<int>(first_edge.from->id());
    int ref_last_pent = static_cast<int>(second_edge.from->id());


    auto candidates = find_L_reductions(G, ref_size, ref_first_pent, ref_last_pent);
    if (candidates.empty()) {
        return true;
    }

    std::uint32_t ref_x2 = edge_neighborhood_code(first_edge, use_next);
    std::vector<l_reduction> stage = candidates;

    {
        std::vector<l_reduction> next;
        next.reserve(stage.size());
        for (const auto& r : stage) {
            std::uint32_t v = edge_neighborhood_code(r.first_edge, r.use_next);
            if (v < ref_x2) {
                return false;
            }
            if (v == ref_x2) {
                next.push_back(r);
            }
        }
        stage.swap(next);
    }

    if (stage.empty()) {
        return true;
    }

    std::uint32_t ref_x3 = edge_neighborhood_code(second_edge, use_next);
    {
        std::vector<l_reduction> next;
        next.reserve(stage.size());
        for (const auto& r : stage) {
            std::uint32_t v = edge_neighborhood_code(r.second_edge, r.use_next);
            if (v < ref_x3) {
                return false;
            }
            if (v == ref_x3) {
                next.push_back(r);
            }
        }
        stage.swap(next);
    }

    if (stage.empty()) {
        return true;
    }


    std::uint32_t ref_x4 = path_neighborhood_code(*this, 7);
    {
        std::vector<l_reduction> next;
        next.reserve(stage.size());
        for (const auto& r : stage) {
            std::uint32_t v = path_neighborhood_code(r, 7);
            if (v < ref_x4) {
                return false;
            }
            if (v == ref_x4) {
                next.push_back(r);
            }
        }
        stage.swap(next);
    }

    if (stage.empty()) {
        return true;
    }

    LCandidate ref_cand;
    ref_cand.start = first_edge;
    ref_cand.use_next = use_next;
    ref_cand.i = ref_size;

    LSignatureState ref_state(G, ref_cand);

    std::vector<LCandidate> cand_data(stage.size());
    std::vector<LSignatureState> states;
    states.reserve(stage.size());
    for (std::size_t i = 0; i < stage.size(); ++i) {
        cand_data[i].start = stage[i].first_edge;
        cand_data[i].use_next = stage[i].use_next;
        cand_data[i].i = stage[i].size;
        states.emplace_back(G, cand_data[i]);
    }

    std::vector<bool> alive(states.size(), true);
    std::size_t alive_count = states.size();
    std::size_t prefix_len = 0;

    while (alive_count > 0) {
        if (!ref_state.finished()) {
            ref_state.extend_step();
        }
        for (std::size_t i = 0; i < states.size(); ++i) {
            if (!alive[i]) {
                continue;
            }
            if (!states[i].finished()) {
                states[i].extend_step();
            }
        }

        const auto& ref_sig = ref_state.signature();
        std::size_t ref_len = ref_sig.size();

        bool any_progress = false;

        for (std::size_t i = 0; i < states.size(); ++i) {
            if (!alive[i]) {
                continue;
            }

            const auto& sig = states[i].signature();
            std::size_t cand_len = sig.size();

            std::size_t j = prefix_len;
            if (j >= ref_len && j >= cand_len) {
                continue;
            }

            any_progress = true;

            for (; j < ref_len || j < cand_len; ++j) {
                bool ref_has = j < ref_len;
                bool cand_has = j < cand_len;

                if (!ref_has && !cand_has) {
                    break;
                }

                if (!ref_has && cand_has) {
                    if (ref_state.finished()) {
                        alive[i] = false;
                        --alive_count;
                    }
                    break;
                }

                if (ref_has && !cand_has) {
                    if (states[i].finished()) {
                        return false;
                    }
                    break;
                }

                int vr = ref_sig[j];
                int vc = sig[j];
                if (vc < vr) {
                    return false;
                }
                if (vc > vr) {
                    alive[i] = false;
                    --alive_count;
                    break;
                }
            }
        }

        prefix_len = ref_len;

        if (!any_progress) {
            break;
        }
    }

    return true;
}
