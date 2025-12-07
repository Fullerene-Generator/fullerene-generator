#include "expansions/l_expansion.h"
#include <expansions/l0_signature_state.h>
#include <queue>

void build_L0_rails(const dual_fullerene& G,
    const directed_edge& e0,
    bool use_next,
    std::array<int, 3>& path,
    std::array<int, 3>& para)
{
    auto e = e0;
    for (int k = 0; k < 3; ++k) {
        path[k] = (int)e.from->id();
        auto einv = e.inverse();
        auto side = use_next ? einv.prev_around() : einv.next_around();
        para[k] = (int)side.to()->id();
        e = use_next ? e.right_turn(3) : e.left_turn(3);
    }
}

std::vector<L0Candidate> find_L0_candidates(const dual_fullerene& G)
{
    std::vector<L0Candidate> out;

    for (auto node : G.get_nodes_5()) {
        for (std::size_t i = 0; i < node->degree(); ++i) {
            directed_edge e{ node, i };

            for (bool use_next : { true, false }) {
                std::array<int, 3> P{}, Q{};
                build_L0_rails(G, e, use_next, P, Q);

                if (G.get_node((unsigned)Q[2])->degree() == 5)
                    out.push_back({ e, use_next, P, Q });
            }
        }
    }
    return out;
}


bool L0_expansion::validate() const
{
    return G_.get_node((unsigned)cand_.para[2])->degree() == 5;
}

void L0_expansion::apply() const
{
    const auto& c = cand_;

    const int u0 = c.path[0];
    const int u1 = c.path[1];
    const int u2 = c.path[2];

    const int w0 = c.para[0];
    const int w1 = c.para[1];
    const int w2 = c.para[2];

    const int h1 = G_.add_vertex(node_type::NODE_6);
    const int h2 = G_.add_vertex(node_type::NODE_6);

    G_.move_neighbourhood(u0, h1);
    G_.move_neighbourhood(w2, h2);

    auto u0_node = G_.get_node(u0);
    auto w2_node = G_.get_node(w2);
    auto u1_node = G_.get_node(u1);
    auto u2_node = G_.get_node(u2);
    auto w0_node = G_.get_node(w0);
    auto w1_node = G_.get_node(w1);
    auto h1_node = G_.get_node(h1);
    auto h2_node = G_.get_node(h2);

    if (c.use_next) {
        G_.add_neighbour_after(h1, u1, u0);
        G_.add_neighbour_after(h2, w1, w2);

        u0_node->add_neighbor(u1_node);
        u0_node->add_neighbor(w2_node);
        u0_node->add_neighbor(w1_node);
        u0_node->add_neighbor(w0_node);
        u0_node->add_neighbor(h1_node);

        w2_node->add_neighbor(w1_node);
        w2_node->add_neighbor(u0_node);
        w2_node->add_neighbor(u1_node);
        w2_node->add_neighbor(u2_node);
        w2_node->add_neighbor(h2_node);
    }
    else {
        G_.add_neighbour_before(h1, u1, u0);
        G_.add_neighbour_before(h2, w1, w2);

        u0_node->add_neighbor(h1_node);
        u0_node->add_neighbor(w0_node);
        u0_node->add_neighbor(w1_node);
        u0_node->add_neighbor(w2_node);
        u0_node->add_neighbor(u1_node);

        w2_node->add_neighbor(h2_node);
        w2_node->add_neighbor(u2_node);
        w2_node->add_neighbor(u1_node);
        w2_node->add_neighbor(u0_node);
        w2_node->add_neighbor(w1_node);
    }

    G_.replace_neighbour(u1, w1, w2);
    G_.replace_neighbour(u1, w0, u0);
    G_.replace_neighbour(u2, w1, w2);
    G_.replace_neighbour(w0, u1, u0);
    G_.replace_neighbour(w1, u1, u0);
    G_.replace_neighbour(w1, u2, w2);
}

std::vector<std::unique_ptr<base_expansion>>
find_L0_expansions(dual_fullerene& G)
{
    std::vector<std::unique_ptr<base_expansion>> out;

    const auto candidates = find_L0_candidates(G);
    std::size_t n = candidates.size();
    if (n == 0) {
        return out;
    }

    std::vector<L0SignatureState> states;
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
            auto e = std::make_unique<L0_expansion>(G, candidates[i]);
            if (e->validate()) {
                out.push_back(std::move(e));
            }
        }
    }

    return out;
}

