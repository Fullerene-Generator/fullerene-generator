#include <expansions/l_reduction.h>
#include <iostream>

std::vector<LReduction>
find_L_reductions(const dual_fullerene& G,
    int size,
    int skip_pent_1,
    int skip_pent_2)
{
    std::vector<LReduction> out;
    if (size < 1) {
        return out;
    }

    int path_len = size + 1;
    for (const auto& pent : G.get_nodes_5()) {
        auto start_node = pent;
        unsigned int start_id = start_node->id();
        int deg = static_cast<int>(start_node->degree());

        for (int i = 0; i < deg; ++i) {
            directed_edge e0{ start_node, static_cast<std::size_t>(i) };
            for (bool use_next : { true, false }) {
                int outside_hex1 = use_next ? e0.prev_around(2).to()->degree() : e0.next_around(2).to()->degree();
                
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
                int outside_hex2 = use_next ? e.next_around(1).to()->degree() : e.prev_around(1).to()->degree();
  
                if (outside_hex2 != 6) {
                    continue;
                }

                int first_pent = path[0];
                int last_pent = path[path_len - 1];

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

                LReduction r;
                r.first_edge = e0;
                r.second_edge = second_edge;
                r.use_next = use_next;
                r.size = size;

                out.push_back(std::move(r));
            }
        }
    }

    return out;
}
