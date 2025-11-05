#ifndef EXPANSIONS_H
#define EXPANSIONS_H

#include "fullerene.h"
#include <vector>
#include <array>
using namespace std;

struct L0Candidate {
    DirEdge start;
    bool use_next;
    array<int, 3> path;
    array<int, 3> para;
};

struct L0Result {
    int new_hex_1;
    int new_hex_2;
    array<int, 3> path;
    array<int, 3> para;
    bool use_next;
};

void build_L0_rails(const dual_fullerene& G, DirEdge e0, bool use_next, int path[3], int para[3]);
vector<L0Candidate> find_L0_candidates(const dual_fullerene& G);
void apply_L0(dual_fullerene& G, const L0Candidate& c);

#endif
