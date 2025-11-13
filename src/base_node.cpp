#include <fullerene/base_node.h>

void base_node::add_neighbor(const std::shared_ptr<base_node>& n) {
    neighbors_.push_back(std::weak_ptr(n));
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

std::optional<directed_edge> base_node::get_edge(const std::size_t index = 0) {
    if (index < degree()) {
        return directed_edge(shared_from_this(), index);
    }
    return std::nullopt;
}

std::optional<directed_edge> base_node::get_edge(const std::shared_ptr<const base_node>& neighbor) {
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        if (neighbors_.at(i).lock() == neighbor) {
            return directed_edge(shared_from_this(), i);
        }
    }
    return std::nullopt;
}