#include "expansions.h"

void build_L0_rails(const dual_fullerene& G, half_edge e0, bool use_next, int path[3], int para[3]) {
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
            half_edge e{ u,i };
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

L0Result apply_L0(dual_fullerene& G, const L0Candidate& c) {
    int u0 = c.path[0], u1 = c.path[1], u2 = c.path[2];
    int w0 = c.para[0], w1 = c.para[1], w2 = c.para[2];

    int H = G.add_vertex(node_type::NODE_6);
    int K = G.add_vertex(node_type::NODE_6);

    G.add_edge_between(H, H, u0, w0);

    G.remove_edge_between(u1, w0);
    G.add_edge_between(u1, u0, H, u0);
    G.add_edge_between(H, u1, w0, u0);

    G.add_edge_between(H, w0, w1, w0);

    G.add_edge_between(K, K, w2, u2);

    G.remove_edge_between(u2, w1);
    G.add_edge_between(u2, u1, K, w2);
    G.add_edge_between(K, u2, w1, w2);

    G.remove_edge_between(u1, w1);
    G.add_edge_between(K, w1, u1, H);

    G.add_edge_between(H, w1, K, u1);

    assert(G.degree(u0) == 5 && G.degree(w2) == 5);
    assert(G.degree(H) == 6 && G.degree(K) == 6);

    return { H, K, {u0,u1,u2}, {w0,w1,w2}, c.use_next };
}
