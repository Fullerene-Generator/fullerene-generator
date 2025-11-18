#include <fullerene/base_node.h>

void base_node::resize_(const std::size_t neighbors) {
    neighbors_.resize(neighbors);
    edges_.resize(neighbors);
}

void base_node::set_neighbor_at(const std::size_t index, const std::shared_ptr<base_node>& n) {
    neighbors_[index] = std::weak_ptr(n);
}

std::size_t base_node::degree() const {
    return neighbors_.size();
}

std::shared_ptr<base_node> base_node::neighbor_at(const std::size_t index) const {
    return neighbors_.at(index).lock();
}

const std::vector<std::weak_ptr<base_node>>& base_node::neighbors() const {
    return neighbors_;
}

directed_edge base_node::get_edge(const std::size_t index = 0) {
    if (index < degree()) {
        return { shared_from_this(), index };
    }

    throw std::out_of_range("Node " + std::to_string(id_) + " doesn't have a neighbor with index " + std::to_string(index) );
}

directed_edge base_node::get_edge(const std::shared_ptr<const base_node>& other) {
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        if (neighbors_.at(i).lock() == other) {
            return { shared_from_this(), i };
        }
    }

    throw std::invalid_argument("Node " + std::to_string(other->id()) + " is not a neighbor of node " + std::to_string(id_) );
}

bool base_node::is_neighbor_of(const std::shared_ptr<base_node> &other) const {
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        if (neighbors_.at(i).lock() == other) {
            return true;
        }
    }

    return false;
}

edge_data& base_node::get_edge_data(const std::size_t index) {
    return edges_.at(index);
}

void base_node::clear_all_edge_data() {
    for (int i = 0; i < degree(); i++) {
        edges_[i] = {};
    }
}
