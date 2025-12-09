#include "fullerene/base_node.h"
#include <stdexcept>
#include <string>
#include <algorithm>


void base_node::add_neighbor(const std::shared_ptr<base_node>& n) {
    neighbors_.push_back(std::weak_ptr<base_node>(n));
    edges_.push_back({});
}

void base_node::resize_(const std::size_t neighbors) {
    neighbors_.resize(neighbors);
    edges_.resize(neighbors);
}

void base_node::set_neighbor_at(const std::size_t index, const std::shared_ptr<base_node>& n) {
    neighbors_[index] = std::weak_ptr(n);
}

void base_node::add_neighbour_after(const std::shared_ptr<base_node>& after,
    const std::shared_ptr<base_node>& new_n) 
{
    auto it = std::find_if(neighbors_.begin(), neighbors_.end(),
        [&](const std::weak_ptr<base_node>& w) {
            auto sp = w.lock();
            return sp && sp == after;
        });

    if (it == neighbors_.end())
        throw std::invalid_argument("Node " + std::to_string(after->id()) +
            " is not a neighbor of node " + std::to_string(id_));

    const auto idx = static_cast<std::size_t>(std::distance(neighbors_.begin(), it));
    neighbors_.insert(it + 1, std::weak_ptr<base_node>(new_n));
    edges_.push_back({});
}

void base_node::add_neighbour_before(const std::shared_ptr<base_node>& before,
    const std::shared_ptr<base_node>& new_n) {
    auto it = std::find_if(neighbors_.begin(), neighbors_.end(),
        [&](const std::weak_ptr<base_node>& w) {
            auto sp = w.lock();
            return sp && sp == before;
        });

    if (it == neighbors_.end())
        throw std::invalid_argument("Node " + std::to_string(before->id()) +
            " is not a neighbor of node " + std::to_string(id_));

    const auto idx = static_cast<std::size_t>(std::distance(neighbors_.begin(), it));
    neighbors_.insert(it, std::weak_ptr<base_node>(new_n));
    edges_.push_back({});
}

void base_node::remove_neighbor(const std::shared_ptr<base_node>& n) {
    auto it = std::find_if(neighbors_.begin(), neighbors_.end(),
        [&](const std::weak_ptr<base_node>& w) {
            auto sp = w.lock();
            return sp && sp == n;
        });

    if (it == neighbors_.end())
        throw std::invalid_argument("Node " + std::to_string(n->id()) +
            " is not a neighbor of node " + std::to_string(id_));

    const auto idx = static_cast<std::size_t>(std::distance(neighbors_.begin(), it));
    neighbors_.erase(it);
    edges_.erase(edges_.begin() + static_cast<std::ptrdiff_t>(idx));
}

void base_node::replace_neighbor(const std::shared_ptr<base_node>& old_n,
    const std::shared_ptr<base_node>& new_n) {
    auto it = std::find_if(neighbors_.begin(), neighbors_.end(),
        [&](const std::weak_ptr<base_node>& w) {
            auto sp = w.lock();
            return sp && sp == old_n;
        });

    if (it == neighbors_.end())
        throw std::invalid_argument("Node " + std::to_string(old_n->id()) +
            " is not a neighbor of node " + std::to_string(id_));

    *it = std::weak_ptr<base_node>(new_n);
}

void base_node::move_neighborhood_from(const std::shared_ptr<base_node>& other) {
    neighbors_ = std::move(other->neighbors_);
    edges_ = std::move(other->edges_);
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

directed_edge base_node::get_edge(const std::size_t index) {
    if (index < degree()) {
        return { shared_from_this(), index };
    }

    throw std::out_of_range("Node " + std::to_string(id_) + " doesn't have a neighbor with index " +
        std::to_string(index));
}

directed_edge base_node::get_edge(const std::shared_ptr<const base_node>& other) {
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        if (neighbors_.at(i).lock() == other) {
            return { shared_from_this(), i };
        }
    }

    throw std::invalid_argument("Node " + std::to_string(other->id()) +
        " is not a neighbor of node " + std::to_string(id_));
}

bool base_node::is_neighbor_of(const std::shared_ptr<base_node>& other) const {
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
    for (std::size_t i = 0; i < degree(); i++) {
        edges_[i] = {};
    }
}
