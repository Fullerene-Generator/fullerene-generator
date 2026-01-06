#ifndef L_SIGNATURE_STATE_H
#define L_SIGNATURE_STATE_H

#include <fullerene/dual_fullerene.h>
#include <expansions/l_expansion.h>
#include <vector>

class signature_state {
    const dual_fullerene* graph_;
    const expansion_candidate* candidate_;
    std::vector<int> signature_;
    std::vector<unsigned int> bfs_order_;
    std::vector<directed_edge> base_edges_;
    std::vector<int> index_of_;
    std::size_t bfs_front_;
    bool finished_;
    int color_offset_;

public:
    signature_state(const dual_fullerene& G, const expansion_candidate& c);

    void extend_step();
    [[nodiscard]] bool finished() const;
    [[nodiscard]] const std::vector<int>& signature() const;
};

#endif