#ifndef B_REDUCTION_H
#define B_REDUCTION_H

#include <vector>

#include <expansions/base_reduction.h>
#include <expansions/b_expansion.h>

struct b_reduction final : base_reduction {
    int length_pre_bend = 0;
    int length_post_bend = 0;

    [[nodiscard]] int x0() const override;
    [[nodiscard]] int x1() const override;

    void apply(dual_fullerene& G, const expansion_candidate& c) const override;
};

std::vector<b_reduction>
find_b_reductions(const dual_fullerene& G, int length_pre_bend, int length_post_bend, int skip_pent, int skip_index, bool skip_clockwise);

#endif
