#include "expansions.h"

void build_L0_rails(const dual_fullerene& G, DirEdge e0, bool use_next, int path[3], int para[3]) {
    auto e = e0;
    for (int k = 0; k < 3; ++k) {
        path[k] = e.u;
        auto einv = G.invers(e);
        auto side = use_next ? G.prev_around(einv) : G.next_around(einv);
        para[k] = G.neighbor_id(side.u, side.i);
        e = use_next ? G.straight_step_next(e) : G.straight_step_prev(e);
    }
}

vector<L0Candidate> find_L0_candidates(const dual_fullerene& G) {
    vector<L0Candidate> out;
    for (int u : G.degree5_vertices()) {
        for (uint8_t i = 0; i < 5; ++i) {
            DirEdge e{ u,i };
            for (bool use_next : {true, false}) {
                int P[3], Q[3];
                build_L0_rails(G, e, use_next, P, Q);
                if (G.degree(Q[2]) == 5) {
                    out.push_back({ e,use_next,{P[0],P[1],P[2]},{Q[0],Q[1],Q[2]} });
                }
            }
        }
    }
    return out;
}

void apply_L0(dual_fullerene& G, const L0Candidate& c) {
    int u0 = c.path[0], u1 = c.path[1], u2 = c.path[2];
    int w0 = c.para[0], w1 = c.para[1], w2 = c.para[2];

    int h1 = G.add_vertex(node_type::NODE_6);
    int h2 = G.add_vertex(node_type::NODE_6);
    G.move_neighbourhood(u0, h1);
    G.move_neighbourhood(w2, h2);
    if (c.use_next) {
        G.add_neighbour_after(h1, u1, u0);
        G.add_neighbour_after(h2, w1, w2);
        G.nodes[u0]->neighbors = { G.nodes[u1], G.nodes[w2], G.nodes[w1], G.nodes[w0], G.nodes[h1] };
        G.nodes[w2]->neighbors = { G.nodes[w1], G.nodes[u0], G.nodes[u1], G.nodes[u2], G.nodes[h2] };
    }
    else {
        G.add_neighbour_before(h1, u1, u0);
        G.add_neighbour_before(h2, w1, w2);
        G.nodes[u0]->neighbors = { G.nodes[h1], G.nodes[w0], G.nodes[w1], G.nodes[w2], G.nodes[u1] };
        G.nodes[w2]->neighbors = { G.nodes[h2], G.nodes[u2], G.nodes[u1], G.nodes[u0], G.nodes[w1] };
    }

    G.replace_neighbour(u1, w1, w2);
    G.replace_neighbour(u1, w0, u0);
    G.replace_neighbour(u2, w1, w2);
    G.replace_neighbour(w0, u1, u0);
    G.replace_neighbour(w1, u1, u0);
    G.replace_neighbour(w1, u2, w2);
}
