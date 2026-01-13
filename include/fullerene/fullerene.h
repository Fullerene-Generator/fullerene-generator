#ifndef FULLERENE_H
#define FULLERENE_H
#include <array>
#include <vector>
#include <string>
#include <cmath>

class fullerene {
    std::string id_;
    std::string parent_id_;
    bool is_ipr_;
    std::vector<std::array<unsigned int, 3>> adjacency_;
    std::array<unsigned int, 5> outer_face_nodes_;

public:
    explicit fullerene(const std::string& id,
                        const std::string& parent_id,
                        const bool is_ipr,
                        const std::vector<std::array<unsigned int, 3>>& adjacency,
                        const std::array<unsigned int, 5> &outer_face):
                        id_(id),
                        parent_id_(parent_id),
                        is_ipr_(is_ipr),
                        adjacency_(adjacency),
                        outer_face_nodes_(outer_face) {};

    [[nodiscard]] std::vector<std::array<unsigned int, 3>> get_adjacency() const { return adjacency_; }
    [[nodiscard]] size_t get_size() const { return adjacency_.size(); }
    [[nodiscard]] std::array<unsigned int, 5> get_outer_face_nodes() const { return outer_face_nodes_; }
    [[nodiscard]] std::string write_all() const noexcept;
    std::string get_parent_id() const { return parent_id_; }
    friend std::ostream &operator<<(std::ostream &os, const fullerene &f);
};

#endif //FULLERENE_H
