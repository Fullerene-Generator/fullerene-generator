#ifndef L_REDUCTION_H
#define L_REDUCTION_H

#include <vector>

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <expansions/l_expansion.h>

struct l_reduction {
    directed_edge first_edge;
    directed_edge second_edge;
    bool use_next;
    int size;

    [[nodiscard]] bool is_canonical(const dual_fullerene& G, int min_size) const;
    void apply(dual_fullerene& G, const l_expansion_candidate& c) const;
};

std::vector<l_reduction>
find_l_reductions(const dual_fullerene& G,
    int size,
    int skip_pent_1 = -1,
    int skip_pent_2 = -1);

#endif
