#ifndef L_REDUCTION_H
#define L_REDUCTION_H

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <vector>

struct LReduction {
    directed_edge first_edge;
    directed_edge second_edge;
    bool use_next;
    int size;
};

std::vector<LReduction>
find_L_reductions(const dual_fullerene& G,
    int size,
    int skip_pent_1 = -1,
    int skip_pent_2 = -1);

#endif
