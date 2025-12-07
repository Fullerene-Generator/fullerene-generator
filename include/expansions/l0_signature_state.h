#ifndef L0_SIGNATURE_STATE_H
#define L0_SIGNATURE_STATE_H

#include <fullerene/dual_fullerene.h>
#include <expansions/l_expansion.h>
#include <vector>

class L0SignatureState {
    const dual_fullerene* graph_;
    const L0Candidate* candidate_;
    std::vector<int> signature_;
    std::vector<unsigned int> bfs_order_;
    std::vector<directed_edge> base_edges_;
    std::vector<int> index_of_;
    std::size_t bfs_front_;
    bool finished_;

public:
    L0SignatureState(const dual_fullerene& G, const L0Candidate& c);

    void extend_step();
    bool finished() const;
    const std::vector<int>& signature() const;
};

#endif