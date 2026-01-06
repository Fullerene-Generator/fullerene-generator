#include <queue>
#include <embeddings/embedder.h>
#include <embeddings/forces.h>

std::vector<unsigned int> embedder::compute_bfs_depth(const graph &f) {
    auto depth = std::vector<unsigned int>(f.adjacency.size());
    auto visited = std::vector<bool>(f.adjacency.size(), false);

    std::queue<unsigned> q;

    for (const auto v : f.outer) {
        depth.at(v) = 0;
        visited.at(v) = true;
        q.push(v);
    }

    while (!q.empty()) {
        auto v = q.front();
        q.pop();

        for (const auto u : f.adjacency.at(v)) {
            if (!visited.at(u)) {
                depth.at(u) = depth.at(v) + 1;
                visited.at(u) = true;
                q.push(u);
            }
        }
    }

    return depth;
}

void embedder::find_pentagons_starting_at(std::vector<std::array<unsigned, 5>> &pentagons,
    std::array<unsigned, 5> &current_pentagon, const std::vector<std::array<unsigned, 3>> &adjacency,
    const unsigned starting_node, const unsigned current_node, const unsigned depth, std::vector<bool> &visited) {

    current_pentagon.at(depth) = current_node;

    if (depth == 4) {
        for (const auto n : adjacency[current_node]) {
            if (n == starting_node && current_pentagon.at(1) < current_pentagon.at(4)) {
                pentagons.push_back(current_pentagon);
            }
        }

        return;
    }

    visited.at(current_node) = true;

    for (const auto n : adjacency[current_node]) {
        if (n < starting_node || visited.at(n)) {
            continue;
        }

        find_pentagons_starting_at(pentagons, current_pentagon, adjacency, starting_node, n, depth + 1, visited);
    }

    visited.at(current_node) = false;
}

std::vector<std::array<unsigned, 5>> embedder::find_pentagons(const graph &f) {
    auto& adjacency = f.adjacency;
    const auto n = adjacency.size();

    auto pentagons = std::vector<std::array<unsigned, 5>>();
    auto pentagon = std::array<unsigned, 5>();
    auto visited = std::vector(n, false);

    pentagons.reserve(12);

    for (int i = 0; i < n; i++) {
        find_pentagons_starting_at(pentagons, pentagon, adjacency, i, i, 0, visited);
    }

    return pentagons;
}

std::vector<std::array<double, 2>> embedder::compute_tutte(const graph& f) {
    auto& adjacency = f.adjacency;
    auto& outer = f.outer;
            
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

std::vector<std::array<double,3>> embedder::compute_spectral_realization(const graph& g) {
    const auto n = static_cast<long long>(g.adjacency.size());

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

    for (long long v = 0; v < n; ++v) {
        for (const unsigned u : g.adjacency[v]) {
            A(v, u) = 1.0;
        }
    }

    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

    if (es.info() != Eigen::Success) {
        throw std::runtime_error("Self adjoint solver returned with an error");
    }

    const auto& V = es.eigenvectors().reverse();

    std::vector<std::array<double, 3>> embedding(n);

    for (int v = 0; v < n; ++v) {
        embedding[v] = {
            V(v, 1),
            V(v, 2),
            V(v, 3)
        };
    }

    double max_len = 0;

    for (unsigned int v = 0; v < n; v++) {
        max_len = std::max(max_len,
            sqrt(embedding[v][0] * embedding[v][0]
            + embedding[v][1] * embedding[v][1]
            + embedding[v][2] * embedding[v][2]));
    }

    if (max_len > 0) {
        for (unsigned int v = 0; v < n; v++) {
            for (int i = 0; i < 3; i++) {
                embedding[v][i] /= max_len;
            }
        }
    }

    return embedding;
}

std::vector<std::array<double, 3>> embedder::compute_tutte_sphere_mapping(const graph &f) {
    const auto n = static_cast<long long>(f.adjacency.size());

    auto depth = compute_bfs_depth(f);
    const auto max_depth = *std::ranges::max_element(depth);

    const auto tutte_embedding = compute_tutte(f);
    std::array<double, 2> barycenter = {0, 0};

    for (unsigned int v = 0; v < n; v++) {
        barycenter[0] += tutte_embedding[v][0];
        barycenter[1] += tutte_embedding[v][1];
    }

    barycenter[0] /= static_cast<double>(n);
    barycenter[1] /= static_cast<double>(n);

    std::vector<std::array<double, 3>> embedding(n);

    for (int v = 0; v < n; v++) {
        const auto phi = (depth[v] + 0.5) * M_PI / (max_depth + 1);
        const auto theta = atan2(tutte_embedding[v][1] - barycenter[1], tutte_embedding[v][0] - barycenter[0]);

        embedding[v][0] = std::sin(phi) * std::cos(theta);
        embedding[v][1] = std::sin(phi) * std::sin(theta);
        embedding[v][2] = std::cos(phi);
    }

    return embedding;
}

std::vector<std::array<double, 3>> embedder::compute_3d_force_embedding(const graph &f) {
    auto embedding = compute_tutte_sphere_mapping(f);

    force_params params;
    relax_bond_springs(f, embedding, params);

    return embedding;
}
