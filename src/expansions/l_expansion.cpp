#include "expansions/l_expansion.h"
#include <expansions/l_signature_state.h>
#include <queue>
#include <iostream>
#include <unordered_set>

void build_l_rails(const dual_fullerene& G,
    const directed_edge& start,
    bool clockwise,
    int length,
    std::vector<int>& path,
    std::vector<int>& parallel_path)
{
    const int len = length + 3;
    path.resize(len);
    parallel_path.resize(len);
    auto e = start;
    for (int k = 0; k < len; ++k) {
        path[k] = (int)e.from->id();
        auto e_inverse = e.inverse();
        auto side = clockwise ? e_inverse.prev_around() : e_inverse.next_around();
        parallel_path[k] = (int)side.to()->id();
        e = clockwise ? e.right_turn(3) : e.left_turn(3);
    }
}

static bool l_patch_vertices_unique(const std::vector<int>& path,
    const std::vector<int>& para)
{
    std::unordered_set<int> seen;
    seen.reserve(path.size() + para.size());

    for (int v : path) {
        if (!seen.insert(v).second) {
            return false;
        }
    }
    for (int v : para) {
        if (!seen.insert(v).second) {
            return false;
        }
    }
    return true;
}

std::vector<l_expansion_candidate> find_l_candidates(const dual_fullerene& G, int x)
{
    std::vector<l_expansion_candidate> out;

    for (const auto& node : G.get_nodes_5()) {
        for (int i = 0; i < node->degree(); ++i) {
            directed_edge e{ node, static_cast<std::size_t>(i) };

            for (bool use_next : { true, false }) {
                std::vector<int> P, Q;
                build_l_rails(G, e, use_next, x, P, Q);

                if ((G.get_node((unsigned)Q[x+2])->degree() == 5) && l_patch_vertices_unique(P, Q))
                    out.push_back({ e, use_next, x, std::move(P), std::move(Q) });
            }
        }
    }
    return out;
}


bool l_expansion::validate() const
{
    return G_.get_node((unsigned)cand_.parallel_path[cand_.length + 2])->degree() == 5;
}

void l_expansion::apply()
{
    const auto& c = cand_;
    int i = cand_.length;
    const int u0 = c.path[0];
    const int u1 = c.path[1];
    const int w0 = c.parallel_path[0];
    const int w1 = c.parallel_path[1];

    //first hexagon outside
    const int h1 = G_.add_vertex(node_type::NODE_6);
    auto h1_node = G_.get_node(h1);
    
    G_.move_neighbourhood(u0, h1);

    auto u0_node = G_.get_node(u0);;
    auto u1_node = G_.get_node(u1);
    auto w0_node = G_.get_node(w0);
    auto w1_node = G_.get_node(w1);
    

    if (c.clockwise) {
        G_.add_neighbour_after(h1, u1, u0);
        u0_node->add_neighbor(u1_node);
        u0_node->add_neighbor(w1_node);
        u0_node->add_neighbor(w0_node);
        u0_node->add_neighbor(h1_node);
    }
    else {
        G_.add_neighbour_before(h1, u1, u0);
        u0_node->add_neighbor(h1_node);
        u0_node->add_neighbor(w0_node);
        u0_node->add_neighbor(w1_node);
        u0_node->add_neighbor(u1_node);
    }

    G_.replace_neighbour(u1, w0, u0);
    G_.replace_neighbour(w0, u1, u0);
    G_.replace_neighbour(w1, u1, u0);

    //loop for corridor of hexagons
    int corridor_v = u0;

    for (int j = 0; j < i; j++) {
        int h = G_.add_vertex(node_type::NODE_6);
        int u_first = c.path[j + 1];
        int u_second = c.path[j + 2];
        int w_first = c.parallel_path[j + 1];
        int w_second = c.parallel_path[j + 2];
        auto u_first_node = G_.get_node(u_first);
        auto u_second_node = G_.get_node(u_second);
        auto w_first_node = G_.get_node(w_first);
        auto w_second_node = G_.get_node(w_second);
        auto corridor_node = G_.get_node(corridor_v);
        auto h_node = G_.get_node(h);
        h_node->clear_neighbors();


        if (c.clockwise) {
            G_.add_neighbour_after(corridor_v, u_first, h);
            h_node->add_neighbor(w_first_node);
            h_node->add_neighbor(corridor_node);
            h_node->add_neighbor(u_first_node);
            h_node->add_neighbor(u_second_node);
            h_node->add_neighbor(w_second_node);
        }
        else {
            G_.add_neighbour_before(corridor_v, u_first, h);
            h_node->add_neighbor(w_second_node);
            h_node->add_neighbor(u_second_node);
            h_node->add_neighbor(u_first_node);
            h_node->add_neighbor(corridor_node);
            h_node->add_neighbor(w_first_node);
        }

        if (j == 0) {
            inv_first_ = corridor_node->get_edge(h_node);
        }
        
        G_.replace_neighbour(u_first, w_first, h);
        G_.replace_neighbour(w_first, u_second, h);
        G_.replace_neighbour(u_second, w_first, h);
        G_.replace_neighbour(w_second, u_second, h);
        corridor_v = h;

    }
    //second hexagon outside
    const int h2 = G_.add_vertex(node_type::NODE_6);
    auto h2_node = G_.get_node(h2);
    int u_first = c.path[i + 1];
    int u_second = c.path[i + 2];
    int w_first = c.parallel_path[i + 1];
    int w_second = c.parallel_path[i + 2];
    auto u_first_node = G_.get_node(u_first);
    auto u_second_node = G_.get_node(u_second);
    auto w_first_node = G_.get_node(w_first);
    auto w_second_node = G_.get_node(w_second);
    auto corridor_node = G_.get_node(corridor_v);

    G_.move_neighbourhood(w_second, h2);

    if (c.clockwise) {
        G_.add_neighbour_after(corridor_v, u_first, w_second);
        G_.add_neighbour_after(h2, w_first, w_second);
        w_second_node->add_neighbor(w_first_node);
        w_second_node->add_neighbor(corridor_node);
        w_second_node->add_neighbor(u_first_node);
        w_second_node->add_neighbor(u_second_node);
        w_second_node->add_neighbor(h2_node);
    }
    else {
        G_.add_neighbour_before(corridor_v, u_first, w_second);
        G_.add_neighbour_before(h2, w_first, w_second);
        w_second_node->add_neighbor(h2_node);
        w_second_node->add_neighbor(u_second_node);
        w_second_node->add_neighbor(u_first_node);
        w_second_node->add_neighbor(corridor_node);
        w_second_node->add_neighbor(w_first_node);
    }

    if (i == 0) {
        inv_first_ = corridor_node->get_edge(w_second_node);
    }
    inv_second_ = w_second_node->get_edge(corridor_node);
    G_.replace_neighbour(u_first, w_first, w_second);
    G_.replace_neighbour(u_second, w_first, w_second);
    G_.replace_neighbour(w_first, u_second, w_second);

}

std::vector<std::unique_ptr<base_expansion>>
find_l_expansions(dual_fullerene& G, int i)
{
    std::vector<std::unique_ptr<base_expansion>> out;

    const auto candidates = find_l_candidates(G,i);
    std::size_t n = candidates.size();
    if (n == 0) {
        return out;
    }

    std::vector<l_signature_state> states;
    states.reserve(n);
    for (const auto& c : candidates) {
        states.emplace_back(G, c);
    }

    struct Group {
        std::vector<std::size_t> members;
        std::size_t prefix_len;
    };

    std::queue<Group> groups;
    Group initial;
    initial.members.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        initial.members.push_back(i);
    }
    initial.prefix_len = 0;
    groups.push(initial);

    std::vector<bool> is_representative(n, false);

    while (!groups.empty()) {
        Group g = groups.front();
        groups.pop();

        if (g.members.empty()) {
            continue;
        }

        if (g.members.size() == 1) {
            std::size_t idx = g.members.front();
            if (!is_representative[idx]) {
                is_representative[idx] = true;
            }
            continue;
        }

        for (std::size_t idx : g.members) {
            states[idx].extend_step();
        }

        std::vector<std::vector<std::size_t>> subgroups;
        std::vector<std::size_t> reps;

        for (std::size_t idx : g.members) {
            const auto& sig = states[idx].signature();

            bool placed = false;
            for (std::size_t k = 0; k < reps.size(); ++k) {
                std::size_t rep_idx = reps[k];
                const auto& rep_sig = states[rep_idx].signature();

                if (sig.size() != rep_sig.size()) {
                    continue;
                }

                bool equal = true;
                std::size_t start = g.prefix_len;
                std::size_t end = sig.size();
                for (std::size_t p = start; p < end; ++p) {
                    if (sig[p] != rep_sig[p]) {
                        equal = false;
                        break;
                    }
                }

                if (equal) {
                    subgroups[k].push_back(idx);
                    placed = true;
                    break;
                }
            }

            if (!placed) {
                reps.push_back(idx);
                subgroups.push_back(std::vector<std::size_t>{idx});
            }
        }

        for (std::size_t k = 0; k < subgroups.size(); ++k) {
            auto& members = subgroups[k];
            if (members.empty()) {
                continue;
            }

            if (members.size() == 1) {
                std::size_t idx = members.front();
                if (!is_representative[idx]) {
                    is_representative[idx] = true;
                }
                continue;
            }

            bool all_finished = true;
            for (std::size_t idx : members) {
                if (!states[idx].finished()) {
                    all_finished = false;
                    break;
                }
            }

            if (all_finished) {
                std::size_t idx = members.front();
                if (!is_representative[idx]) {
                    is_representative[idx] = true;
                }
            }
            else {
                Group ng;
                ng.members = std::move(members);
                ng.prefix_len = states[reps[k]].signature().size();
                groups.push(std::move(ng));
            }
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (is_representative[i]) {
            auto e = std::make_unique<l_expansion>(G, candidates[i]);
            if (e->validate()) {
                out.push_back(std::move(e));
            }
        }
    }

    return out;
}
