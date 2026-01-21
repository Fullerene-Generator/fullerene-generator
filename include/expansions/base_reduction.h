#ifndef BASE_REDUCTION_H
#define BASE_REDUCTION_H

#include <cstdint>
#include <memory>
#include <vector>

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <expansions/base_expansion.h>

class base_reduction {
public:
    directed_edge first_edge;
    directed_edge second_edge;
    bool use_next;

    virtual ~base_reduction() = default;

    [[nodiscard]] virtual int x0() const = 0;
    [[nodiscard]] virtual int x1() const { return -x0(); }

    [[nodiscard]] bool is_canonical(const dual_fullerene& G, int min_x0, int l1, int l2) const;

    virtual void apply(dual_fullerene& G, const expansion_candidate& c) const = 0;

protected:
    [[nodiscard]] std::uint32_t x2_code() const;
    [[nodiscard]] std::uint32_t x3_code() const;
    [[nodiscard]] std::uint32_t x4_code() const;

    void fill_signature_candidate(expansion_candidate& out) const;
};

std::vector<std::unique_ptr<base_reduction>>
find_all_reductions(const dual_fullerene& G, int x0, int skip_pent, int skip_index, bool skip_clockwise, int skip_l1, int skip_l2);

int limit_by_reduction_distances(const dual_fullerene& G, int cur_best);


#endif
