#include <expansions/b_expansion.h>

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
        auto e_inverse = e.inverse();
        auto side = clockwise ? e_inverse.prev_around() : e_inverse.next_around();
        parallel_path[k - 1] = static_cast<int>(side.to()->id());
        e = clockwise ? e.right_turn(3) : e.left_turn(3);
    }

    e = clockwise ? e.next_around() : e.prev_around();

    for (int k = length_pre_bend + 2; k < total_length + 1; ++k) {
        path[k] = static_cast<int>(e.from->id());
        auto e_inverse = e.inverse();
        auto side = clockwise ? e_inverse.prev_around() : e_inverse.next_around();
        parallel_path[k - 1] = static_cast<int>(side.to()->id());
        e = clockwise ? e.right_turn(3) : e.left_turn(3);
    }

    path[total_length + 1] = static_cast<int>(e.to()->id());
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

                if ((G.get_node(static_cast<unsigned>(Q[Q.size() - 1]))->degree() == 5) && patch_nodes_unique(P, Q))
                    out.push_back({ e, clockwise, length_pre_bend, length_post_bend, std::move(P), std::move(Q) });
            }
        }
    }

    return out;
}
