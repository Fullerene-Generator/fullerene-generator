#include "fullerene.h"

// ---- base_node implementation ----

void base_node::add_neighbor(const std::shared_ptr<base_node>& n) {
    neighbors_.push_back(std::weak_ptr(n));
}

std::size_t base_node::neighbor_count() const {
    return neighbors_.size();
}

// ---- dual_fullerene implementation ----

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
            throw std::invalid_argument("Pentagon node " + std::to_string(i) + " has degree " +
                std::to_string(deg) + ", expected 5");
        if (i >= 12 && deg != 6)
            throw std::invalid_argument("Hexagon node " + std::to_string(i) + " has degree " +
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


