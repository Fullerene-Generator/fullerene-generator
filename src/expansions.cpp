#include "fullerene/expansions.h"

void build_L0_rails(const dual_fullerene& G,
    const directed_edge& e0,
    const bool use_next,
    std::array<int, 3>& path,
    std::array<int, 3>& para)
{
    auto e = e0;

    for (int k = 0; k < 3; ++k) {
        path[k] = static_cast<int>(e.from->id());

        auto einv = e.inverse();
        auto side = use_next ? einv.prev_around() : einv.next_around();

        para[k] = static_cast<int>(side.to()->id());
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
                std::array<int, 3> P{};
                std::array<int, 3> Q{};

                build_L0_rails(G, e, use_next, P, Q);

                if (G.get_node(static_cast<unsigned int>(Q[2]))->degree() == 5) {
                    out.push_back(L0Candidate{ e, use_next, P, Q });
                }
            }
        }
    }

    return out;
}

void apply_L0(dual_fullerene& G, const L0Candidate& c)
{
    const int u0 = c.path[0];
    const int u1 = c.path[1];
    const int u2 = c.path[2];

    const int w0 = c.para[0];
    const int w1 = c.para[1];
    const int w2 = c.para[2];

    const int h1 = G.add_vertex(node_type::NODE_6);
    const int h2 = G.add_vertex(node_type::NODE_6);

    G.move_neighbourhood(u0, h1);
    G.move_neighbourhood(w2, h2);

    auto u0_node = G.get_node(static_cast<unsigned int>(u0));
    auto w2_node = G.get_node(static_cast<unsigned int>(w2));

    auto u1_node = G.get_node(static_cast<unsigned int>(u1));
    auto u2_node = G.get_node(static_cast<unsigned int>(u2));
    auto w0_node = G.get_node(static_cast<unsigned int>(w0));
    auto w1_node = G.get_node(static_cast<unsigned int>(w1));
    auto h1_node = G.get_node(static_cast<unsigned int>(h1));
    auto h2_node = G.get_node(static_cast<unsigned int>(h2));

    if (c.use_next) {
        G.add_neighbour_after(h1, u1, u0);
        G.add_neighbour_after(h2, w1, w2);

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
        G.add_neighbour_before(h1, u1, u0);
        G.add_neighbour_before(h2, w1, w2);

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

    G.replace_neighbour(u1, w1, w2);
    G.replace_neighbour(u1, w0, u0);
    G.replace_neighbour(u2, w1, w2);
    G.replace_neighbour(w0, u1, u0);
    G.replace_neighbour(w1, u1, u0);
    G.replace_neighbour(w1, u2, w2);
}
