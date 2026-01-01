#ifndef L_EXPANSION_H
#define L_EXPANSION_H

#include <fullerene/dual_fullerene.h>
#include <fullerene/directed_edge.h>
#include <expansions/base_expansion.h>
#include <utility>
#include <vector>

struct l_candidate {
    directed_edge start;
    bool use_next;
    int i;
    std::vector<int> path;
    std::vector<int> para;
};

void build_l_rails(const dual_fullerene& G,
    const directed_edge& e0,
    bool use_next,
    int i,
    std::vector<int>& path,
    std::vector<int>& para);

std::vector<l_candidate> find_l_candidates(const dual_fullerene& G, int i);

class l_expansion final : public base_expansion {
    l_candidate cand_;
    directed_edge inv_first_;
    directed_edge inv_second_;

public:
    explicit l_expansion(dual_fullerene& G, l_candidate  c)
        : base_expansion(G), cand_(std::move(c)) {
    }

    [[nodiscard]] bool validate() const override;
    void apply() override;

    [[nodiscard]] directed_edge inverse_first_edge() const { return inv_first_; }
    [[nodiscard]] directed_edge inverse_second_edge() const { return inv_second_; }
    [[nodiscard]] const l_candidate& candidate() const { return cand_; }
};

std::vector<std::unique_ptr<base_expansion>>
find_l_expansions(dual_fullerene& G, int i);

#endif
