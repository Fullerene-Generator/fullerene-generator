#ifndef FULLERENE_GENERATOR_EMBEDDER_H
#define FULLERENE_GENERATOR_EMBEDDER_H
#include <set>
#include <Eigen/Dense>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif //M_PI

using angle_key = std::tuple<double, double, double>;

struct graph {
    std::vector<std::array<unsigned,3>> adjacency;
    std::array<unsigned,5> outer;
};

class embedder {
    static std::vector<unsigned> compute_bfs_depth(const graph& f);
    static void find_pentagons_starting_at(std::vector<std::array<unsigned, 5>>& pentagons,
        std::array<unsigned, 5>& current_pentagon,const std::vector<std::array<unsigned,3>>& adjacency,
        unsigned starting_node, unsigned current_node, unsigned depth, std::vector<bool>& visited);
    static std::vector<std::array<unsigned, 5>> find_pentagons(const graph& f);
    static std::set<angle_key> find_pentagon_angles(const graph& f);

public:
    inline static angle_key make_angle_key (double i, double j, double k);
    static std::vector<std::array<double, 2>> compute_tutte(const graph& f);
    static std::vector<std::array<double, 3>> compute_spectral_realization(const graph& f);
    static std::vector<std::array<double, 3>> compute_tutte_sphere_mapping(const graph& f);
    static std::vector<std::array<double, 3>> compute_3d_force_embedding(const graph& f);
};


#endif //FULLERENE_GENERATOR_EMBEDDER_H

