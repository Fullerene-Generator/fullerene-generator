#include <expansions/b_reduction.h>

#include <algorithm>
#include <stdexcept>

std::vector<b_reduction>
find_b_reductions(const dual_fullerene& G, int length_pre_bend, int length_post_bend)
{
    std::vector<b_reduction> out;

    for (const auto& node : G.get_nodes_5()) {
        for (int i = 0; i < node->degree(); ++i) {
            directed_edge e{ node, static_cast<std::size_t>(i) };

            for (bool clockwise : { true, false }) {
                std::vector<int> P, Q;
                build_b_rails(G, e, clockwise, length_pre_bend, length_post_bend, P, Q);

                if (P.empty()) continue;

                if (G.get_node(static_cast<unsigned>(P.back()))->degree() != 5) continue;
                if (!patch_nodes_unique(P, Q)) continue;

                auto last_node = G.get_node(static_cast<unsigned>(P.back()));
                auto prev_node = G.get_node(static_cast<unsigned>(P[P.size() - 2]));
                directed_edge second_edge = last_node->get_edge(prev_node);

                b_reduction r;
                r.first_edge = e;
                r.second_edge = second_edge;
                r.use_next = clockwise;
                r.length_pre_bend = length_pre_bend;
                r.length_post_bend = length_post_bend;

                out.push_back(r);
            }
        }
    }

    return out;
}

int b_reduction::x0() const
{
    return length_pre_bend + length_post_bend + 2;
}

int b_reduction::x1() const
{
    return -(std::max(length_pre_bend, length_post_bend) + 1);
}



void b_reduction::apply(dual_fullerene& G, const expansion_candidate& c_base) const
{
    const auto& c = static_cast<const b_expansion_candidate&>(c_base);

    const int i1 = c.length_pre_bend;
    const int i2 = c.length_post_bend;
    const int i_total = i1 + i2;

    const int created = i_total + 3;
    const int h1 = static_cast<int>(G.total_nodes()) - created;
    int h = static_cast<int>(G.total_nodes()) - 1;

    // undo second external hexagon 
    {
        const int h2 = h;

        const int u_first = c.path[i_total + 3];
        const int u_second = c.path[i_total + 4];
        const int w_first = c.parallel_path[i_total + 1];
        const int w_second = c.parallel_path[i_total + 2];

        auto h2_node = G.get_node(static_cast<unsigned>(h2));
        auto u_second_node = G.get_node(static_cast<unsigned>(u_second));

        h2_node->remove_neighbor(u_second_node);
        G.move_neighborhood(h2, u_second);
        G.pop_last_node6();
        --h;

        G.replace_neighbor(u_first, u_second, w_second);
        G.replace_neighbor(w_first, u_second, u_first);
        G.replace_neighbor(w_second, u_second, u_first);
    }

    // undo corridor post-bend
    for (int j = i_total; j >= i1 + 1; --j) {
        const int u_first = c.path[j + 2];
        const int u_second = c.path[j + 3];
        const int w_first = c.parallel_path[j];
        const int w_second = c.parallel_path[j + 1];

        G.replace_neighbor(w_first, h, u_first);
        G.replace_neighbor(u_first, h, w_second);
        G.replace_neighbor(w_second, h, u_first);
        G.replace_neighbor(u_second, h, w_second);

        G.pop_last_node6();
        --h;
    }

    // undo bend hexagon
    {
        const int u_first = c.path[i1 + 1];
        const int u_second = c.path[i1 + 2];
        const int u_third = c.path[i1 + 3];
        const int w_first = c.parallel_path[i1 + 1];

        G.replace_neighbor(u_first, h, w_first);
        G.replace_neighbor(w_first, h, u_second);
        G.replace_neighbor(u_second, h, w_first);
        G.replace_neighbor(u_third, h, w_first);

        G.pop_last_node6();
        --h;
    }

    // undo corridor pre-bend
    for (int j = i1; j >= 1; --j) {
        const int u_first = c.path[j];
        const int u_second = c.path[j + 1];
        const int w_first = c.parallel_path[j];
        const int w_second = c.parallel_path[j + 1];

        G.replace_neighbor(u_first, h, w_first);
        G.replace_neighbor(w_first, h, u_second);
        G.replace_neighbor(u_second, h, w_first);
        G.replace_neighbor(w_second, h, u_second);

        G.pop_last_node6();
        --h;
    }

    // undo first external hexagon
    {
        const int u0 = c.path[0];
        const int u1 = c.path[1];
        const int w0 = c.parallel_path[0];
        const int w1 = c.parallel_path[1];

        auto h1_node = G.get_node(static_cast<unsigned>(h1));
        auto u0_node = G.get_node(static_cast<unsigned>(u0));

        h1_node->remove_neighbor(u0_node);
        G.move_neighborhood(h1, u0);
        G.pop_last_node6();

        G.replace_neighbor(u1, u0, w0);
        G.replace_neighbor(w0, u0, u1);
        G.replace_neighbor(w1, u0, u1);
    }
}
