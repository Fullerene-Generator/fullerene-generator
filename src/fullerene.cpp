#include <fstream>
#include <ranges>
#include <fullerene/fullerene.h>
#include <Eigen/Dense>



void fullerene::compute_tutte_embedding() {
    const auto n = static_cast<long long>(adjacency_.size());
    embedding_2d_.resize(n, { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() });

    std::array<std::array<double, 2>,5> outer_face_coords{};

    for (int i = 0; i < 5; i++) {
        const double angle = 2.0 * M_PI * i / 5.0 + 0.5 * M_PI;
        outer_face_coords[i] = { std::cos(angle), std::sin(angle) };
    }

    for (int i = 0; i < 5; i++) {
        const unsigned int v = outer_face_nodes_[i];
        embedding_2d_[v] = outer_face_coords[i];
    }

    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(n);

    for (unsigned int v = 0; v < n; v++) {
        const auto &neighbors = adjacency_[v];

        if (embedding_2d_[v][0] != std::numeric_limits<double>::max()) {
            L(v,v) = 1.0;
            bx(v) = embedding_2d_[v][0];
            by(v) = embedding_2d_[v][1];
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
        embedding_2d_[v][0] = x(v);
        embedding_2d_[v][1] = y(v);
    }
}

std::string fullerene::write_graph() const noexcept {
    std::ostringstream ss;

    // fullerene size
    ss << adjacency_.size() << "\n";

    // nodes of the outer face
    for (unsigned v : outer_face_nodes_) {
        ss << v << ' ';
    }
    ss << "\n";

    // adjacency of vertices
    for (auto const& adj : adjacency_) {
        ss << adj[0] << " " << adj[1] << " " << adj[2] << "\n";
    }

    return ss.str();
}

std::string fullerene::write_embedding() const noexcept {
    std::stringstream ss;

    // vertex coordinates
    for (const auto& coords : embedding_2d_) {
        ss << coords[0] << ' ' << coords[1] << '\n';
    }

    return ss.str();
}

std::string fullerene::write_all() const noexcept {
    std::stringstream ss;

    ss << write_graph();
    if (has_2d_embedding()) ss << write_embedding();

    return ss.str();
}

fullerene fullerene::read_graph(std::istream& is) {
    std::size_t n;
    is >> n;

    std::array<unsigned,5> outer{};
    for (int i = 0; i < 5; i++) {
        is >> outer[i];
    }

    std::vector<std::array<unsigned,3>> adj(n);
    for (std::size_t i = 0; i < n; i++) {
        is >> adj[i][0] >> adj[i][1] >> adj[i][2];
    }

    return fullerene(adj, outer);
}

std::ostream & operator<<(std::ostream &os, const fullerene &f) {
    os << f.write_all();
    return os;
}
