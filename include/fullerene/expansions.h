#ifndef EXPANSIONS_H
#define EXPANSIONS_H

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <array>
#include <vector>

struct L0Candidate {
    directed_edge             start;
    bool                      use_next;
    std::array<int, 3>        path;
    std::array<int, 3>        para;
};

void build_L0_rails(const dual_fullerene& G,
    const directed_edge& e0,
    bool use_next,
    std::array<int, 3>& path,
    std::array<int, 3>& para);

std::vector<L0Candidate> find_L0_candidates(const dual_fullerene& G);
void apply_L0(dual_fullerene& G, const L0Candidate& c);

#endif
