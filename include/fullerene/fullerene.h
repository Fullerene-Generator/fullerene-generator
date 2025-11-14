#ifndef FULLERENE_H
#define FULLERENE_H
#include <array>
#include <vector>

class fullerene {
    std::array<unsigned int, 5> outer_face_nodes_;
    std::vector<std::array<unsigned int, 3>> adjacency_;

public:
    explicit fullerene(const std::vector<std::array<unsigned int, 3>>& adjacency,
                       const std::array<unsigned int, 5> outer_face): adjacency_(adjacency),
                       outer_face_nodes_(outer_face) {}
};

#endif //FULLERENE_H
