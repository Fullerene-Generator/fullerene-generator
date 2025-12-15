#include"expansions/f_expansion.h"

bool f_expansion::validate() const {
    for (const auto& u : v_->neighbors()) {
        if (u.lock()->type() == node_type::NODE_6) {
            return false;
        }
    }

    return true;
}

void f_expansion::apply(){
    std::array<std::shared_ptr<node_6>, 5> new_nodes{};

    for (auto & new_node : new_nodes) {
        new_node = node_6::create_sized(G_.total_nodes());
        G_.add_node(new_node);
    }

    for (int i = 0; i < 5; i++) {
        auto& v1 = new_nodes[i];
        auto& v2 = new_nodes[(i + 1) % 5];

        v1->set_neighbor_at(3, v2);
        v2->set_neighbor_at(0, v1);
    }

    std::array<directed_edge, 5> edges;

    for (int i = 0; i < 5; i++) {
        edges[i] = v_->get_edge(i);
    }

    for (int i = 0; i < 5; i++) {
        auto edge_l = edges[i].left_turn(2);
        auto edge_r = edges[i].right_turn(2);
        auto edge_l_inv = edge_l.inverse();
        auto edge_r_inv = edge_r.inverse();

        const auto& curr_new_node = new_nodes[i];
        const auto& next_new_node = new_nodes[(i + 1) % 5];

        edge_l.change_destination(curr_new_node);
        curr_new_node->set_neighbor_at(4, edge_l.from);

        edge_l_inv.change_destination(next_new_node);
        next_new_node->set_neighbor_at(1, edge_l_inv.from);

        edge_r.change_destination(next_new_node);
        next_new_node->set_neighbor_at(5, edge_r.from);

        edge_r_inv.change_destination(next_new_node);
        next_new_node->set_neighbor_at(2, edge_r_inv.from);
    }
}