#ifndef L_EXPANSION_H
#define L_EXPANSION_H

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <expansions/base_expansion.h>
#include <vector>


struct LCandidate {
    directed_edge start;
    bool use_next;
    int i;
    std::vector<int> path;
    std::vector<int> para;
};

void build_L_rails(const dual_fullerene& G,
    const directed_edge& e0,
    bool use_next,
    int i,
    std::vector<int>& path,
    std::vector<int>& para);

std::vector<LCandidate> find_L_candidates(const dual_fullerene& G, int i);


class L_expansion final : public base_expansion {
    LCandidate cand_;

public:
    explicit L_expansion(dual_fullerene& G, const LCandidate& c)
        : base_expansion(G), cand_(c) {
    }

    [[nodiscard]] bool validate() const override;
    void apply() const override;

    const LCandidate& candidate() const { return cand_; }
};

std::vector<std::unique_ptr<base_expansion>>
find_L_expansions(dual_fullerene& G, int i);

#endif
