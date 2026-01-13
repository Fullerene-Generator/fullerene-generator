#include <queue>
#include <expansions/b_expansion.h>
#include <expansions/signature_state.h>

void build_b_rails(const dual_fullerene& G,
                   const directed_edge& start,
                   bool clockwise,
                   int length_pre_bend,
                   int length_post_bend,
                   std::vector<int>& path,
                   std::vector<int>& parallel_path)
{
    const int total_length = length_pre_bend + length_post_bend + 3;
    path.resize(total_length + 2);
    parallel_path.resize(total_length);

    auto e = start;
    path[0] = static_cast<int>(e.from->id());
    e = clockwise ? e.right_turn(3) : e.left_turn(3);

    for (int k = 1; k < length_pre_bend + 2; ++k) {
        path[k] = static_cast<int>(e.from->id());
        auto side = clockwise ? e.next_around(2) : e.prev_around(2);
        parallel_path[k - 1] = static_cast<int>(side.to()->id());
        e = clockwise ? e.right_turn(3) : e.left_turn(3);
    }

    e = clockwise ? e.next_around() : e.prev_around();

    for (int k = length_pre_bend + 2; k < total_length + 1; ++k) {
        path[k] = static_cast<int>(e.from->id());
        auto side = clockwise ? e.next_around() : e.prev_around();
        parallel_path[k - 1] = static_cast<int>(side.to()->id());
        e = clockwise ? e.right_turn(3) : e.left_turn(3);
    }

    path[total_length + 1] = static_cast<int>(e.from->id());
}

std::vector<b_expansion_candidate> find_b_candidates(const dual_fullerene& G,
    int length_pre_bend,
    int length_post_bend)
{
    std::vector<b_expansion_candidate> out;

    for (const auto& node : G.get_nodes_5()) {
        for (int i = 0; i < node->degree(); ++i) {
            directed_edge e{ node, static_cast<std::size_t>(i) };

            for (bool clockwise : { true, false }) {
                std::vector<int> P, Q;
                build_b_rails(G, e, clockwise, length_pre_bend, length_post_bend, P, Q);

                if ((G.get_node(static_cast<unsigned>(P[P.size() - 1]))->degree() == 5) && patch_nodes_unique(P, Q))
                    out.push_back({ e, clockwise, std::move(P), std::move(Q), length_pre_bend, length_post_bend });
            }
        }
    }

    return out;
}

bool b_expansion::validate() const {
    return G_.get_node(static_cast<unsigned>(cand_.path[cand_.path.size() - 1]))->degree() == 5;
}

void b_expansion::apply() {
    const auto& c = cand_;
    int i1 = cand_.length_pre_bend;
    int i2 = cand_.length_post_bend;
    int i_total = i1 + i2;
    const int u0 = c.path[0];
    const int u1 = c.path[1];
    const int w0 = c.parallel_path[0];
    const int w1 = c.parallel_path[1];

    // first external hexagon
    const int h1 = G_.add_vertex(node_type::NODE_6);
    auto h1_node = G_.get_node(h1);

    G_.move_neighborhood(u0, h1);

    auto u0_node = G_.get_node(u0);;
    auto u1_node = G_.get_node(u1);
    auto w0_node = G_.get_node(w0);
    auto w1_node = G_.get_node(w1);

    if (c.clockwise) {
        G_.add_neighbor_after(h1, u1, u0);
        u0_node->add_neighbor(u1_node);
        u0_node->add_neighbor(w1_node);
        u0_node->add_neighbor(w0_node);
        u0_node->add_neighbor(h1_node);
    }
    else {
        G_.add_neighbor_before(h1, u1, u0);
        u0_node->add_neighbor(h1_node);
        u0_node->add_neighbor(w0_node);
        u0_node->add_neighbor(w1_node);
        u0_node->add_neighbor(u1_node);
    }

    G_.replace_neighbor(u1, w0, u0);
    G_.replace_neighbor(w0, u1, u0);
    G_.replace_neighbor(w1, u1, u0);

    // corridor of hexagons (pre-bend)
    int corridor_v = u0;

    {
        for (int j = 1; j <= i1; j++) {
            int h = G_.add_vertex(node_type::NODE_6);
            int u_first = c.path[j];
            int u_second = c.path[j + 1];
            int w_first = c.parallel_path[j];
            int w_second = c.parallel_path[j + 1];
            auto u_first_node = G_.get_node(u_first);
            auto u_second_node = G_.get_node(u_second);
            auto w_first_node = G_.get_node(w_first);
            auto w_second_node = G_.get_node(w_second);
            auto corridor_node = G_.get_node(corridor_v);
            auto h_node = G_.get_node(h);
            h_node->clear_neighbors();

            if (c.clockwise) {
                G_.add_neighbor_after(corridor_v, u_first, h);
                h_node->add_neighbor(w_first_node);
                h_node->add_neighbor(corridor_node);
                h_node->add_neighbor(u_first_node);
                h_node->add_neighbor(u_second_node);
                h_node->add_neighbor(w_second_node);
            }
            else {
                G_.add_neighbor_before(corridor_v, u_first, h);
                h_node->add_neighbor(w_second_node);
                h_node->add_neighbor(u_second_node);
                h_node->add_neighbor(u_first_node);
                h_node->add_neighbor(corridor_node);
                h_node->add_neighbor(w_first_node);
            }

            G_.replace_neighbor(u_first, w_first, h);
            G_.replace_neighbor(w_first, u_second, h);
            G_.replace_neighbor(u_second, w_first, h);
            G_.replace_neighbor(w_second, u_second, h);
            corridor_v = h;
        }
    }

    // bend hexagon
    {
        int h = G_.add_vertex(node_type::NODE_6);
        int u_first = c.path[i1 + 1];
        int u_second = c.path[i1 + 2];
        int u_third = c.path[i1 + 3];
        int w_first = c.parallel_path[i1 + 1];
        auto u_first_node = G_.get_node(u_first);
        auto u_second_node = G_.get_node(u_second);
        auto u_third_node = G_.get_node(u_third);
        auto w_first_node = G_.get_node(w_first);
        auto corridor_node = G_.get_node(corridor_v);
        auto h_node = G_.get_node(h);
        h_node->clear_neighbors();

        if (c.clockwise) {
            G_.add_neighbor_after(corridor_v, u_first, h);
            h_node->add_neighbor(w_first_node);
            h_node->add_neighbor(corridor_node);
            h_node->add_neighbor(u_first_node);
            h_node->add_neighbor(u_second_node);
            h_node->add_neighbor(u_third_node);
        }
        else {
            G_.add_neighbor_before(corridor_v, u_first, h);
            h_node->add_neighbor(u_third_node);
            h_node->add_neighbor(u_second_node);
            h_node->add_neighbor(u_first_node);
            h_node->add_neighbor(corridor_node);
            h_node->add_neighbor(w_first_node);
        }

        G_.replace_neighbor(u_first, w_first, h);
        G_.replace_neighbor(w_first, u_second, h);
        G_.replace_neighbor(u_second, w_first, h);
        G_.replace_neighbor(u_third, w_first, h);
        corridor_v = h;
    }

    // corridor of hexagons (post-bend)
    for (int j = i1 + 1; j <= i_total; j++) {
        int h = G_.add_vertex(node_type::NODE_6);
        int u_first = c.path[j + 2];
        int u_second = c.path[j + 3];
        int w_first = c.parallel_path[j];
        int w_second = c.parallel_path[j + 1];
        auto u_first_node = G_.get_node(u_first);
        auto u_second_node = G_.get_node(u_second);
        auto w_first_node = G_.get_node(w_first);
        auto w_second_node = G_.get_node(w_second);
        auto corridor_node = G_.get_node(corridor_v);
        auto h_node = G_.get_node(h);
        h_node->clear_neighbors();

        if (c.clockwise) {
            G_.add_neighbor_after(corridor_v, u_first, h);
            h_node->add_neighbor(w_first_node);
            h_node->add_neighbor(corridor_node);
            h_node->add_neighbor(u_first_node);
            h_node->add_neighbor(u_second_node);
            h_node->add_neighbor(w_second_node);
        }
        else {
            G_.add_neighbor_before(corridor_v, u_first, h);
            h_node->add_neighbor(w_second_node);
            h_node->add_neighbor(u_second_node);
            h_node->add_neighbor(u_first_node);
            h_node->add_neighbor(corridor_node);
            h_node->add_neighbor(w_first_node);
        }

        G_.replace_neighbor(w_first, u_first, h);
        G_.replace_neighbor(u_first, w_second, h);
        G_.replace_neighbor(w_second, u_first, h);
        G_.replace_neighbor(u_second, w_second, h);
        corridor_v = h;
    }

    // second external hexagon
    {
        const int h2 = G_.add_vertex(node_type::NODE_6);
        auto h2_node = G_.get_node(h2);
        int u_first = c.path[i_total + 3];
        int u_second = c.path[i_total + 4];
        int w_first = c.parallel_path[i_total + 1];
        int w_second = c.parallel_path[i_total + 2];
        auto u_first_node = G_.get_node(u_first);
        auto u_second_node = G_.get_node(u_second);
        auto w_first_node = G_.get_node(w_first);
        auto w_second_node = G_.get_node(w_second);
        auto corridor_node = G_.get_node(corridor_v);

        G_.move_neighborhood(u_second, h2);

        if (c.clockwise) {
            G_.add_neighbor_after(corridor_v, u_first, u_second);
            G_.add_neighbor_after(h2, w_second, u_second);
            u_second_node->add_neighbor(w_first_node);
            u_second_node->add_neighbor(corridor_node);
            u_second_node->add_neighbor(u_first_node);
            u_second_node->add_neighbor(h2_node);
            u_second_node->add_neighbor(w_second_node);
        }
        else {
            G_.add_neighbor_before(corridor_v, u_first, u_second);
            G_.add_neighbor_before(h2, w_second, u_second);
            u_second_node->add_neighbor(h2_node);
            u_second_node->add_neighbor(u_first_node);
            u_second_node->add_neighbor(corridor_node);
            u_second_node->add_neighbor(w_first_node);
            u_second_node->add_neighbor(w_second_node);
        }

        G_.replace_neighbor(u_first, w_second, u_second);
        G_.replace_neighbor(w_first, u_first, u_second);
        G_.replace_neighbor(w_second, u_first, u_second);

        inv_first_ = c.clockwise ? u0_node->get_edge(u1_node).next_around() : u0_node->get_edge(u1_node).prev_around();
        inv_second_ = u_second_node->get_edge(corridor_node);
    }
}

std::vector<std::unique_ptr<base_expansion>> find_b_expansions(dual_fullerene& G,
    int length_pre_bend,
    int length_post_bend)
{
    std::vector<std::unique_ptr<base_expansion>> out;

    const auto candidates = find_b_candidates(G, length_pre_bend, length_post_bend);
    std::size_t n = candidates.size();
    if (n == 0) {
        return out;
    }

    std::vector<signature_state> states;
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
            auto e = std::make_unique<b_expansion>(G, candidates[i]);
            if (e->validate()) {
                out.push_back(std::move(e));
            }
        }
    }

    return out;
}
