#include <fullerene/fullerene.h>
#include <Eigen/Dense>

#include <utility>
#include <algorithm>

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

std::vector<int> dual_fullerene::degree5_vertices() const {
    std::vector<int> out;
    out.reserve(12);
    for (int v = 0; v < (int)nodes.size(); ++v)
        if (degree(v) == 5) out.push_back(v);
    return out;
}

void dual_fullerene::assign_sequential_ids() {
    for (int v = 0; v < (int)nodes.size(); ++v) nodes[v]->id = v;
}

int dual_fullerene::add_vertex(node_type type) {
    int id = static_cast<int>(nodes.size());
    auto v = std::make_shared<node>();
    v->id = id;
    v->type = type;
    v->neighbors.clear();
    nodes.push_back(std::move(v));
    return id;
}

DirEdge dual_fullerene::invers(DirEdge e) const {
    int u = e.u, v = neighbor_id(u, e.i);
    const auto& ring = nodes[v]->neighbors;
    for (uint8_t j = 0; j < ring.size(); ++j)
        if (ring[j].get() == nodes[u].get()) return { v, j };
}

void dual_fullerene::add_neighbour(int v, int nei) {
    nodes[v]->neighbors.push_back(nodes[nei]);
}

void dual_fullerene::add_neighbour_after(int v, int after, int v2) {
    auto it = find(nodes[v]->neighbors.begin(), nodes[v]->neighbors.end(), nodes[after]);
    nodes[v]->neighbors.insert(it + 1, nodes[v2]);
}

void dual_fullerene::add_neighbour_before(int v, int before, int v2) {
    auto it = find(nodes[v]->neighbors.begin(), nodes[v]->neighbors.end(), nodes[before]);
    nodes[v]->neighbors.insert(it, nodes[v2]);
}

void dual_fullerene::remove_edge(int v1, int v2) {
    auto it1 = find(nodes[v1]->neighbors.begin(), nodes[v1]->neighbors.end(), nodes[v2]);
    nodes[v1]->neighbors.erase(it1);
    auto it2 = find(nodes[v2]->neighbors.begin(), nodes[v2]->neighbors.end(), nodes[v1]);
    nodes[v1]->neighbors.erase(it2);
}

void dual_fullerene::replace_neighbour(int v, int old_n, int new_n) {
    auto it = find(nodes[v]->neighbors.begin(), nodes[v]->neighbors.end(), nodes[old_n]);
    *it = nodes[new_n];
}

void dual_fullerene::move_neighbourhood(int from, int to) {
    nodes[to]->neighbors = nodes[from]->neighbors;
    nodes[from]->neighbors.clear();
    for (const auto& x : nodes[to]->neighbors) {
        replace_neighbour(x->id, from, to);
    }
}
