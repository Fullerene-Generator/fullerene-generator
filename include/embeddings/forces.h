#ifndef FULLERENE_GENERATOR_FORCES_H
#define FULLERENE_GENERATOR_FORCES_H
#include <embeddings/embedder.h>
#include <vector>

struct force_params {
    int max_iterations = 5000;
    double step = 0.1;

    double bond_k = 0.5;
    double angle_k = 0.1;
    double radial_k = 0.03;

    double target_bond_length = 1.0;
    double target_angle_hex = 2.0 * M_PI / 3.0;
    double target_angle_pent = 3.0 * M_PI / 5.0;

    double convergence_eps = 1e-6;
};

template <size_t D>
std::array<double, D> barycenter(const std::vector<std::array<double, D>>& pos) {
    std::array<double, D> c{};
    for (const auto& p : pos)
        for (size_t k = 0; k < D; ++k)
            c[k] += p[k];

    for (size_t k = 0; k < D; ++k)
        c[k] /= pos.size();

    return c;
}

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
void apply_angular_forces(const graph &g, std::vector<std::array<double, D>> &pos, std::vector<std::array<double, D>>& force, std::set<angle_key>& pentagon_angles, const force_params &params) {
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
                const double delta = theta - (pentagon_angles.contains(embedder::make_angle_key(i, j, k)) ?
                    params.target_angle_pent :
                    params.target_angle_hex);

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
void apply_radial_repulsion (std::vector<std::array<double, D>> &pos, std::vector<std::array<double, D>>& force, const force_params &params) {
    auto center = barycenter(pos);

    for (size_t i = 0; i < pos.size(); ++i) {
        std::array<double, D> d{};
        double r2 = 0;

        for (size_t k = 0; k < D; ++k) {
            d[k] = pos[i][k] - center[k];
            r2 += d[k] * d[k];
        }

        const double r = std::sqrt(r2);

        if (r == 0.0)
            continue;

        double mag = params.radial_k / (r * r);

        for (size_t k = 0; k < D; ++k)
            force[i][k] += mag * d[k];
    }
}

template <size_t D>
void relax_bond_springs(const graph &g, std::vector<std::array<double, D>> &pos, std::set<angle_key>& pentagon_angles, force_params &params) {
    const auto n = pos.size();

    auto force = std::vector<std::array<double, D>>(n);
    params.target_bond_length = mean_edge_length(g, pos);

    for (int it = 0; it < params.max_iterations; ++it) {
        for (auto& f : force) f.fill(0.0);

        apply_bond_forces(g, pos, force, params);
        apply_angular_forces(g, pos, force, pentagon_angles, params);
        apply_radial_repulsion(pos, force, params);

        auto max_delta = 0.0;

        for (size_t i = 0; i < pos.size(); ++i) {
            auto delta = 0.0;

            for (size_t d = 0; d < D; ++d) {
                auto delta_dim = params.step * force[i][d];
                pos[i][d] += delta_dim;
                delta += delta_dim * delta_dim;
            }

            delta = std::sqrt(delta);
            max_delta = std::max(max_delta, delta);
        }

        if (max_delta <= params.convergence_eps)
            return;
    }
}

#endif //FULLERENE_GENERATOR_FORCES_H
