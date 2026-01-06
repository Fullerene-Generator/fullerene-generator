#ifndef FULLERENE_GENERATOR_FORCES_H
#define FULLERENE_GENERATOR_FORCES_H
#include <embeddings/embedder.h>
#include <vector>

struct force_params {
    int iterations = 1000;
    double step = 0.01;
    double stiffness = 1.0;
};

template <size_t D>
double mean_edge_length(const graph& g, const std::vector<std::array<double, D>>& pos) {
    auto sum = 0.0;
    size_t count = 0;

    for (std::size_t i = 0; i < g.adjacency.size(); ++i) {
        for (auto j : g.adjacency[i]) {
            if (j > i) {
                auto d2 = 0.0;
                for (std::size_t k = 0; k < D; ++k) {
                    const auto diff = pos[i][k] - pos[j][k];
                    d2 += diff * diff;
                }
                sum += std::sqrt(d2);
                count++;
            }
        }
    }
    return sum / static_cast<double>(count);
}

template <size_t D>
void relax_bond_springs(const graph &g, std::vector<std::array<double, D>> &pos, const force_params &params) {
    const auto n = pos.size();
    const double L = mean_edge_length(g, pos);

    auto force = std::vector<std::array<double, D>>(n);

    for (int it = 0; it < params.iterations; ++it) {
        for (auto f: force) {
            f.fill(0);
        }

        for (std::size_t i = 0; i < n; ++i) {
            for (unsigned j: g.adjacency[i]) {
                if (j <= i)
                    continue;

                auto len = 0.0;
                std::array<double, D> d;

                for (std::size_t k = 0; k < D; ++k) {
                    const auto diff = pos[i][k] - pos[j][k];
                    len += diff * diff;
                }

                len = std::sqrt(len);
                if (len == 0)
                    continue;

                auto mag = params.stiffness * (len - L) / len;

                for (std::size_t k = 0; k < D; ++k) {
                    double f = mag * d[k];
                    force[i][k] += f;
                    force[j][k] -= f;
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < D; ++k) {
                pos[i][k] += params.step * force[i][k];
            }
        }
    }
}

#endif //FULLERENE_GENERATOR_FORCES_H