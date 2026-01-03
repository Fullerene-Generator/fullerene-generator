#ifndef FULLERENE_GENERATOR_B_EXPANSION_H
#define FULLERENE_GENERATOR_B_EXPANSION_H
#include <expansions/base_expansion.h>
#include <fullerene/directed_edge.h>

#include <utility>

struct b_expansion_candidate {
    directed_edge start;
    bool clockwise;
    int length_pre_bend;
    int length_post_bend;
    std::vector<int> path;
    std::vector<int> parallel_path;
};

void build_b_rails(const dual_fullerene& G,
    const directed_edge& e0,
    bool clockwise,
    int length_pre_bend,
    int length_post_bend,
    std::vector<int>& path,
    std::vector<int>& parallel_path);

std::vector<b_expansion_candidate> find_b_candidates(const dual_fullerene& G,
    int total_length);

class b_expansion : public base_expansion {
    b_expansion_candidate cand_;

public:
    explicit b_expansion(dual_fullerene &G, b_expansion_candidate c)
        : base_expansion(G), cand_(std::move(c)) { }

    [[nodiscard]] bool validate() const override;
    void apply() override;
};


#endif //FULLERENE_GENERATOR_B_EXPANSION_H