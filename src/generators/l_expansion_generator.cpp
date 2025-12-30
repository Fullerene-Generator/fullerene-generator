#include <generators/l_expansion_generator.h>

namespace {

    l_reduction matching_reduction_from_expansion(const l_expansion& e)
    {
        l_reduction r;
        r.first_edge = e.inverse_first_edge();
        r.second_edge = e.inverse_second_edge();
        r.use_next = e.candidate().use_next;
        r.size = e.candidate().i + 1;
        return r;
    }

} // namespace

void l_expansion_generator::generate(std::size_t up_to)
{
    if (up_to < 20) {
        return;
    }

    {
        auto G = create_c20_fullerene();
        register_and_emit(G);
        dfs_(G, up_to, 3, 1);
        
    }
    
    if (up_to < 28) {
        return;
    }

    {
        auto G = create_c28_fullerene();
        register_and_emit(G);
        dfs_(G, up_to, 4, 1);
    }
    

}



int l_expansion_generator::bound_by_vertex_count_(const dual_fullerene& G, std::size_t up_to)
{
    int primal_v = 20 + 2 * static_cast<int>(G.get_nodes_6().size());
    int dif = static_cast<int>(up_to) - primal_v;
    return dif / 2 - 2;
}



void l_expansion_generator::dfs_(dual_fullerene& G,
    std::size_t up_to,
    int max_i,
    int min_reduction_size)
{
    if (max_i < 0) {
        return;
    }

    std::vector<std::unique_ptr<base_expansion>> expansions;
    expansions.reserve(256);

    for (int i = 0; i <= max_i; ++i) {
        auto v = find_l_expansions(G, i);
        expansions.reserve(expansions.size() + v.size());
        for (auto& up : v) {
            expansions.push_back(std::move(up));
        }
    }

    for (auto& up : expansions) {
        auto* le = static_cast<l_expansion*>(up.get());
        if (!le->validate()) {
            continue;
        }

        const auto cand = le->candidate();

        le->apply();

        const auto red = matching_reduction_from_expansion(*le);

        if (red.is_canonical(G, min_reduction_size)) {
            register_and_emit(G);
            int max_expansion_size = bound_by_vertex_count_(G, up_to);
            int max_i = std::min(max_expansion_size, 8);
            dfs_(G, up_to, max_i, min_reduction_size);
        }

        red.apply(G, cand);
    }
}
