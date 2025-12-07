#include <expansions/l0_signature_state.h>

L0SignatureState::L0SignatureState(const dual_fullerene& G, const L0Candidate& c)
    : graph_(&G),
    candidate_(&c),
    bfs_front_(0),
    finished_(false)
{
    std::size_t n = G.total_nodes();
    index_of_.assign(n, -1);

    signature_.clear();
    bfs_order_.clear();
    base_edges_.clear();

    signature_.reserve(4*n);
    bfs_order_.reserve(n);
    base_edges_.reserve(3*n);


    unsigned int from_id = c.start.from->id();
    unsigned int to_id = c.start.to()->id();

    bfs_order_.push_back(from_id);
    base_edges_.push_back(c.start);
    index_of_[from_id] = 0;

    bfs_order_.push_back(to_id);
    base_edges_.push_back(c.start.inverse());
    index_of_[to_id] = 1;

    signature_.push_back(0);
    signature_.push_back(1);
}

void L0SignatureState::extend_step() {
    if (finished_) {
        return;
    }

    if (bfs_front_ >= bfs_order_.size()) {
        finished_ = true;
        return;
    }

    unsigned int v_id = bfs_order_[bfs_front_];
    directed_edge base_edge = base_edges_[bfs_front_];
    ++bfs_front_;

    auto v_node = graph_->get_node(v_id);
    std::size_t deg = v_node->degree();
    signature_.push_back(static_cast<int>(deg));

    directed_edge e = base_edge;
    for (std::size_t k = 0; k < deg; ++k) {
        auto to_node = e.to();
        unsigned int nid = to_node->id();

        int idx = index_of_[nid];
        if (idx == -1) {
            int new_idx = static_cast<int>(bfs_order_.size());
            index_of_[nid] = new_idx;
            bfs_order_.push_back(nid);
            base_edges_.push_back(e.inverse());
            signature_.push_back(new_idx);
        }
        else {
            signature_.push_back(idx);
        }

        e = e.next_around();
    }

    if (bfs_front_ >= bfs_order_.size()) {
        finished_ = true;
    }
}

bool L0SignatureState::finished() const {
    return finished_;
}

const std::vector<int>& L0SignatureState::signature() const {
    return signature_;
}
