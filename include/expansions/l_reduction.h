#ifndef L_REDUCTION_H
#define L_REDUCTION_H

#include <vector>

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <expansions/base_reduction.h>

struct l_reduction final : base_reduction {
    int size;

    [[nodiscard]] int x0() const override { return size; }
    [[nodiscard]] int x1() const override { return -x0(); }

    void apply(dual_fullerene& G, const expansion_candidate& c) const override;
};

std::vector<l_reduction>
find_l_reductions(const dual_fullerene& G, int size);

#endif
