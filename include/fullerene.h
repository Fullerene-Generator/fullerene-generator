#ifndef FULLERENE_H
#define FULLERENE_H
#include <memory>
#include <vector>

enum class node_type {
    NODE_5,
    NODE_6,
  };

struct node {
    int id;
    node_type type;
    std::vector<std::shared_ptr<node>> neighbors;
};

struct half_edge {
    int      u;    // start vertex id
    uint8_t  i;    // position (slot) in u's clockwise ring: 0..degree(u)-1
};

class dual_fullerene {
    std::vector<std::shared_ptr<node>> nodes;

public:
    explicit dual_fullerene(std::vector<std::shared_ptr<node>> nodes);

    int size() const { return static_cast<int>(nodes.size()); }
    int degree(int u) const { return static_cast<int>(nodes[u]->neighbors.size()); }

    // Id of the neighbor at slot i in u's clockwise ring
    int neighbor_id(int u, uint8_t i) const {
        return nodes[u]->neighbors[i]->id;
    }

    std::vector<int> degree5_vertices() const {
        std::vector<int> out;
        out.reserve(12);
        for (int v = 0; v < (int)nodes.size(); ++v)
            if (degree(v) == 5) out.push_back(v);
        return out;
    }

    void assign_sequential_ids() {
        for (int v = 0; v < (int)nodes.size(); ++v) nodes[v]->id = v;
    }

    // Make a half-edge (u -> neighbors[u][i])
    static half_edge make_half_edge(int u, uint8_t i) { return { u, i }; }

    // NEXT/PREV: move to next/previous slot around the same start vertex (clockwise / counter-clockwise)
    half_edge next_around(half_edge e) const {
        return { e.u, static_cast<uint8_t>((e.i + 1) % degree(e.u)) };
    }
    half_edge prev_around(half_edge e) const {
        return { e.u, static_cast<uint8_t>((e.i + degree(e.u) - 1) % degree(e.u)) };
    }

    half_edge invers(half_edge e) const {
        int u = e.u, v = neighbor_id(u, e.i);
        const auto& ring = nodes[v]->neighbors;
        for (uint8_t j = 0; j < ring.size(); ++j)
            if (ring[j].get() == nodes[u].get()) return { v, j };
    }

    // One straight hop across a triangular face (L0 uses three of these in a row)
    half_edge straight_step_next(half_edge e) const { return next_around(next_around(next_around(invers(e)))); }
    half_edge straight_step_prev(half_edge e) const { return prev_around(prev_around(prev_around(invers(e)))); }


    void replace_neighbor(int u, int oldN, int newN);

    void insert_neighbor_after(int u, int afterN, int newN);

    void add_edge_between(int x, int x_after, int y, int y_after);

    int add_vertex(node_type t);
};

#endif //FULLERENE_H
