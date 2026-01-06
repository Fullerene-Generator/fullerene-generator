#ifndef FULLERENE_GENERATOR_FORCES_H
#define FULLERENE_GENERATOR_FORCES_H
#include <embeddings/embedder.h>
#include <vector>

struct force_params {
    int iterations = 1000;
    double step = 0.01;
    double bond_k = 0.5;
    double angle_k = 0.1;
    double target_bond_length = 1.0;
    double target_angle = 2.0 * M_PI / 3.0;
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
void apply_angular_forces(const graph &g, std::vector<std::array<double, D>> &pos, std::vector<std::array<double, D>>& force, const force_params &params) {
    for (size_t i = 0; i < g.adjacency.size(); ++i) {
        const auto neighbors = g.adjacency[i];

        for (int a = 0; a < 3; ++a) {
            for (int b = a + 1; b < 3; ++b) {
                auto j = neighbors[a];
                auto k = neighbors[b];

                std::array<double, D> u{}, v{};
                double lu = 0, lv = 0;

                for (size_t d = 0; d < D; ++d) {
                    u[d] = pos[j][d] - pos[i][d];
                    v[d] = pos[k][d] - pos[i][d];
                    lu += u[d] * u[d];
                    lv += v[d] * v[d];
                }

                lu = std::sqrt(lu);
                lv = std::sqrt(lv);
                if (lu == 0.0 || lv == 0.0)
                    continue;

                double dot = 0;
                for (size_t d = 0; d < D; ++d)
                    dot += u[d] * v[d];

                double cos_theta = dot / (lu * lv);
                cos_theta = std::clamp(cos_theta, -1.0, 1.0);

                const double theta = std::acos(cos_theta);
                const double delta = theta - params.target_angle;

                double mag = params.angle_k * delta;

                for (size_t d = 0; d < D; ++d) {
                    double fu = mag * (v[d] / lv - cos_theta * u[d] / lu);
                    double fv = mag * (u[d] / lu - cos_theta * v[d] / lv);

                    force[j][d] += fu;
                    force[k][d] += fv;
                    force[i][d] -= (fu + fv);
                }
            }
        }
    }
}

template <size_t D>
void apply_bond_forces(const graph &g, std::vector<std::array<double, D>> &pos, std::vector<std::array<double, D>>& force, const force_params &params) {
    for (std::size_t i = 0; i < g.adjacency.size(); ++i) {
        for (unsigned j: g.adjacency[i]) {
            if (j <= i)
                continue;

            auto len = 0.0;
            std::array<double, D> d;

            for (std::size_t k = 0; k < D; ++k) {
                d[k] = pos[j][k] - pos[i][k];
                len += d[k] * d[k];
            }

            len = std::sqrt(len);
            if (len == 0)
                continue;

            auto mag = params.bond_k * (len - params.target_bond_length) / len;

            for (std::size_t k = 0; k < D; ++k) {
                double f = mag * d[k];
                force[i][k] += f;
                force[j][k] -= f;
            }
        }
    }
}

template <size_t D>
void relax_bond_springs(const graph &g, std::vector<std::array<double, D>> &pos, force_params &params) {
    const auto n = pos.size();

    auto force = std::vector<std::array<double, D>>(n);
    params.target_bond_length = mean_edge_length(g, pos);

    for (int it = 0; it < params.iterations; ++it) {
        for (auto& f : force) f.fill(0.0);

        apply_bond_forces(g, pos, force, params);
        apply_angular_forces(g, pos, force, params);

        for (size_t i = 0; i < pos.size(); ++i)
            for (size_t d = 0; d < D; ++d)
                pos[i][d] += params.step * force[i][d];
    }
}

#endif //FULLERENE_GENERATOR_FORCES_H