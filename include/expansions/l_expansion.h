#ifndef L_EXPANSION_H
#define L_EXPANSION_H

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <expansions/base_expansion.h>
#include <array>
#include <vector>


struct L0Candidate {
    directed_edge start;
    bool use_next;
    std::array<int, 3> path;
    std::array<int, 3> para;
};

void build_L0_rails(const dual_fullerene& G,
    const directed_edge& e0,
    bool use_next,
    std::array<int, 3>& path,
    std::array<int, 3>& para);

std::vector<L0Candidate> find_L0_candidates(const dual_fullerene& G);


class L0_expansion final : public base_expansion {
    L0Candidate cand_;

public:
    explicit L0_expansion(dual_fullerene& G, const L0Candidate& c)
        : base_expansion(G), cand_(c) {
    }

    [[nodiscard]] bool validate() const override;
    void apply() const override;

    const L0Candidate& candidate() const { return cand_; }
};

std::vector<std::unique_ptr<base_expansion>>
find_L0_expansions(dual_fullerene& G);

#endif
