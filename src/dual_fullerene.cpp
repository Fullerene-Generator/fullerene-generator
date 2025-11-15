#include <fullerene/dual_fullerene.h>

template<typename F>
void dual_fullerene::for_each_node(F &&f) const {
    for (const auto& node : nodes_5) f(node);
    for (const auto& node : nodes_6) f(node);
}

dual_fullerene::dual_fullerene(const std::vector<std::vector<unsigned int>>& adjacency) {
    const std::size_t n = adjacency.size();

    std::vector<std::shared_ptr<base_node>> index_to_node(n);

    for (std::size_t i = 0; i < n; i++) {
        if (i < 12) {
            auto p = node_5::create(static_cast<unsigned int>(i));
            index_to_node[i] = p;
            nodes_5.push_back(p);
        } else {
            auto p = node_6::create(static_cast<unsigned int>(i));
            index_to_node[i] = p;
            nodes_6.push_back(p);
        }
    }

    for (std::size_t i = 0; i < n; i++) {
        const auto& neighs = adjacency[i];

        for (const auto j : neighs) {
            if (j >= n)
                throw std::out_of_range("Adjacency index " + std::to_string(j) + " out of range");
            if (j == i)
                throw std::invalid_argument("Self-loop at node " + std::to_string(i));
        }

        const std::size_t deg = neighs.size();
        if (i < 12 && deg != 5)
            throw std::invalid_argument("5-node " + std::to_string(i) + " has degree " +
                std::to_string(deg) + ", expected 5");
        if (i >= 12 && deg != 6)
            throw std::invalid_argument("6-node " + std::to_string(i) + " has degree " +
                std::to_string(deg) + ", expected 6");
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (auto j : adjacency[i]) {
            const auto& neighs = adjacency[j];
            if (std::ranges::find(neighs.begin(), neighs.end(), static_cast<unsigned int>(i)) == neighs.end()) {
                throw std::invalid_argument("Adjacency not symmetric between " +
                    std::to_string(i) + " and " + std::to_string(j));
            }
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (const auto j : adjacency[i]) {
            const auto& a = index_to_node[i];
            const auto& b = index_to_node[j];
            a->add_neighbor(b);
        }
    }
}

fullerene dual_fullerene::to_primal() const {
    const std::size_t V = total_nodes();
    const std::size_t E = (5 * 12 + 6 * (V - 12)) / 2;
    const std::size_t F = E - V + 2;

    std::vector<std::array<unsigned int, 3>> adjacency(F);
    std::vector<unsigned int> counts(F, 0);
    unsigned int face = 0;

    for_each_node([&](const std::shared_ptr<base_node>& node) {
        for (int i = 0; i < node->degree(); i++) {
            auto edge = node->get_edge(i);
            if (edge.data().marked) continue;

            const auto start_node = edge.from;

            do {
                edge.data().marked = true;
                edge.data().rhs_face_index = face;
                edge = edge.right_turn();
            } while (edge.from != start_node);

            face++;
        }
    });

    for_each_node([&](const std::shared_ptr<base_node>& node) {
        for (int i = 0; i < node->degree(); i++) {
            auto edge = node->get_edge(i);

            const auto u = edge.data().rhs_face_index;
            const auto v = edge.inverse().data().rhs_face_index;

            adjacency[u][counts[u]++] = v;
        }
    });

    auto outer_face = nodes_5[0];
    auto outer_face_nodes = std::array<unsigned int, 5>();

    for (int i = 0; i < outer_face->degree(); i++) {
        auto edge = outer_face->get_edge(i).value();
        outer_face_nodes[i] = edge.data().rhs_face_index;
    }

    clear_all_edge_data();

    return std::move(fullerene(adjacency, outer_face_nodes));
}

void dual_fullerene::clear_all_edge_data() const {
    for_each_node([](const std::shared_ptr<base_node>& node) {
        node->clear_all_edge_data();
    });
}
