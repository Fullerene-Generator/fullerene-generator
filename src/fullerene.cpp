#include "fullerene.h"

#include <utility>
#include <algorithm>

dual_fullerene::dual_fullerene(std::vector<std::shared_ptr<node>> nodes) {
    this->nodes = std::move(nodes);
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
