#ifndef EXPANSIONS_H
#define EXPANSIONS_H

#include "fullerene.h"
#include <vector>
#include <array>
using namespace std;

struct L0Candidate {
    half_edge start;
    bool use_next;
    array<int,3> path;
    array<int,3> para;
};

struct L0Result {
    int new_hex_1;
    int new_hex_2;
    array<int, 3> path;
    array<int, 3> para;
    bool use_next;
};

void build_L0_rails(const dual_fullerene& G, half_edge e0, bool use_next, array<int,3>& path, array<int,3>& para)
{
    half_edge e = e0;
    for (int i = 0; i < 3; i++) {
        path[i] = e.u;
        half_edge inv = G.invers(e);
        half_edge side = use_next ? G.prev_around(inv) : G.next_around(inv);
        int c = G.neighbor_id(inv.u, inv.i);
        para[i] = c;
        e = use_next ? G.straight_step_prev(e) : G.straight_step_next(e);
    }
}
vector<L0Candidate> find_L0_candidates(const dual_fullerene& G)
{
    vector<L0Candidate> res;
    for (int u : G.degree5_vertices()) {
        for (int i = 0; i < 5; i++) {
            half_edge e(u, i);
            for (bool use_next : {true, false}) {
                array<int,3> path, para;
                build_L0_rails(G, e, use_next, path, para);
                if (G.degree(para[2]) == 5) {
                    L0Candidate cand = { e, use_next, path, para };
                    res.push_back(cand);
                }
            }
        }
    }
}
L0Result apply_L0(dual_fullerene& G, const L0Candidate& c);

#endif
