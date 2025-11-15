#ifndef FULLERENE_H
#define FULLERENE_H
#include <array>
#include <vector>

class fullerene {
    std::vector<std::array<unsigned int, 3>> adjacency_;

public:
    explicit fullerene(const std::vector<std::array<unsigned int, 3>>& adjacency): adjacency_(adjacency) {}
};

#endif //FULLERENE_H
