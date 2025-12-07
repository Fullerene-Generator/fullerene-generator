#include <embeddings/embedder.h>

std::vector<std::array<double, 2>> embedder::compute_tutte(graph& g) {
    auto& adjacency = g.adjacency;
    auto& outer = g.outer;
            
    const auto n = static_cast<long long>(adjacency.size());
    auto embedding = std::vector<std::array<double, 2>>(n, {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()});

    std::array<std::array<double, 2>,5> outer_face_coords{};

    for (int i = 0; i < 5; i++) {
        const double angle = 2.0 * M_PI * i / 5.0 + 0.5 * M_PI;
        outer_face_coords[i] = { std::cos(angle), std::sin(angle) };
    }

    for (int i = 0; i < 5; i++) {
        const unsigned int v = outer[i];
        embedding[v] = outer_face_coords[i];
    }

    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(n);

    for (unsigned int v = 0; v < n; v++) {
        const auto &neighbors = adjacency[v];

        if (embedding[v][0] != std::numeric_limits<double>::max()) {
            L(v,v) = 1.0;
            bx(v) = embedding[v][0];
            by(v) = embedding[v][1];
            continue;
        }

        L(v,v) = 3.0;

        for (const auto u : neighbors) {
            L(v,u) = -1.0;
        }
    }

    Eigen::VectorXd x = L.colPivHouseholderQr().solve(bx);
    Eigen::VectorXd y = L.colPivHouseholderQr().solve(by);

    for (unsigned int v = 0; v < n; v++) {
        embedding[v][0] = x(v);
        embedding[v][1] = y(v);
    }

    return embedding;
}
