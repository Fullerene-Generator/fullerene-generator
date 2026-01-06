#ifndef FULLERENE_GENERATOR_FORCES_H
#define FULLERENE_GENERATOR_FORCES_H
#include <embeddings/embedder.h>
#include <vector>

enum class force_model {
    none,
    bond_springs
};

struct force_params {
    int iterations = 1000;
    double step = 0.01;
    double stiffness = 1.0;
};

template <size_t D>
double mean_edge_length(const graph& g, const std::vector<std::array<double, D>>& pos);

template <size_t D>
void relax_bond_springs(const graph& g, const std::vector<std::array<double, D>>& pos, const force_params& force);

#endif //FULLERENE_GENERATOR_FORCES_H