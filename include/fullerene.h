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

struct DirEdge {
    int      u;
    uint8_t  i;
};

class dual_fullerene {

public:
    std::vector<std::shared_ptr<node>> nodes;
    explicit dual_fullerene(std::vector<std::shared_ptr<node>> nodes);

    int size() const { return static_cast<int>(nodes.size()); }
    int degree(int u) const { return static_cast<int>(nodes[u]->neighbors.size()); }

    int neighbor_id(int u, uint8_t i) const {
        return nodes[u]->neighbors[i]->id;
    }

    std::vector<int> degree5_vertices() const;

    void assign_sequential_ids();

    int add_vertex(node_type type);

    static DirEdge make_half_edge(int u, uint8_t i) { return { u, i }; }

    DirEdge next_around(DirEdge e) const {
        return { e.u, static_cast<uint8_t>((e.i + 1) % degree(e.u)) };
    }
    DirEdge prev_around(DirEdge e) const {
        return { e.u, static_cast<uint8_t>((e.i + degree(e.u) - 1) % degree(e.u)) };
    }

    DirEdge invers(DirEdge e) const;

    DirEdge straight_step_next(DirEdge e) const { return next_around(next_around(next_around(invers(e)))); }
    DirEdge straight_step_prev(DirEdge e) const { return prev_around(prev_around(prev_around(invers(e)))); }

    void add_neighbour(int v, int nei);
    void add_neighbour_after(int v, int after, int v2);
    void add_neighbour_before(int v, int before, int v2);
    void remove_edge(int v1, int v2);
    void replace_neighbour(int v, int old_n, int new_n);
    void move_neighbourhood(int from, int to);
};

#endif // FULLERENE_H
