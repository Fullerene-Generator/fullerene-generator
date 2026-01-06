#ifndef FULLERENE_GENERATOR_EMBEDDER_H
#define FULLERENE_GENERATOR_EMBEDDER_H
#include <Eigen/Dense>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif //M_PI

struct graph {
    std::vector<std::array<unsigned,3>> adjacency;
    std::array<unsigned,5> outer;
};

class embedder {
    static std::vector<unsigned> compute_bfs_depth(const graph& f);

public:
    static std::vector<std::array<double, 2>> compute_tutte(const graph& f);
    static std::vector<std::array<double, 3>> compute_spectral_realization(const graph& f);
    static std::vector<std::array<double, 3>> compute_tutte_sphere_mapping(const graph& f);
    static std::vector<std::array<double, 3>> compute_3d_force_embedding(const graph& f);
};


#endif //FULLERENE_GENERATOR_EMBEDDER_H

