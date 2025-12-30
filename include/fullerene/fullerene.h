#ifndef FULLERENE_H
#define FULLERENE_H
#include <array>
#include <vector>
#include <string>
#include <cmath>

class fullerene {
    std::vector<std::array<unsigned int, 3>> adjacency_;
    std::array<unsigned int, 5> outer_face_nodes_;

public:
    explicit fullerene(const std::vector<std::array<unsigned int, 3>>& adjacency,
                        const std::array<unsigned int, 5> &outer_face):
                        adjacency_(adjacency),
                        outer_face_nodes_(outer_face) {};

    [[nodiscard]] std::vector<std::array<unsigned int, 3>> get_adjacency() const { return adjacency_; }
    [[nodiscard]] std::array<unsigned int, 5> get_outer_face_nodes() const { return outer_face_nodes_; }
    [[nodiscard]] std::string write_all() const noexcept;
    friend std::ostream &operator<<(std::ostream &os, const fullerene &f);
};

#endif //FULLERENE_H
