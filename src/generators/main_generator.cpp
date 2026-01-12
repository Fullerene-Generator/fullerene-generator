#include <generators/main_generator.h>

#include <expansions/b_expansion.h>
#include <expansions/b_reduction.h>
#include <expansions/l_expansion.h>
#include <expansions/l_reduction.h>
#include <fullerene/construct.h>

#include <algorithm>
#include <memory>
#include <vector>

namespace {

    l_reduction matching_reduction_from_expansion(const l_expansion& e)
    {
        l_reduction r;
        r.first_edge = e.inverse_first_edge();
        r.second_edge = e.inverse_second_edge();
        r.use_next = e.candidate().clockwise;
        r.size = e.candidate().length + 1;
        return r;
    }

    b_reduction matching_reduction_from_expansion(const b_expansion& e)
    {
        b_reduction r;
        r.first_edge = e.inverse_first_edge();
        r.second_edge = e.inverse_second_edge();
        r.use_next = e.candidate().clockwise;
        r.length_pre_bend = e.candidate().length_pre_bend;
        r.length_post_bend = e.candidate().length_post_bend;
        return r;
    }

} 

void main_generator::generate(std::size_t up_to)
{
    if (up_to < 20) {
        return;
    }

    {
        auto G = create_c20_fullerene();
        register_and_emit(G);
        dfs_(G, up_to, 4, 4, 1);
    }

    if (up_to < 28) {
        return;
    }

    {
        auto G = create_c28_fullerene();
        register_and_emit(G);
    }

    for (const auto& [v, c] : counts_) {
        std::cout << v << " vertices: " << c << "\n";
    }
}

int main_generator::bound_by_vertex_count_l(const dual_fullerene& G, std::size_t up_to)
{
    int primal_v = 20 + 2 * static_cast<int>(G.get_nodes_6().size());
    int dif = static_cast<int>(up_to) - primal_v;
    return dif / 2 - 2;
}

int main_generator::bound_by_vertex_count_b(const dual_fullerene& G, std::size_t up_to)
{
    int primal_v = 20 + 2 * static_cast<int>(G.get_nodes_6().size());
    int dif = static_cast<int>(up_to) - primal_v;
    return dif / 2 - 3;
}

void main_generator::dfs_(dual_fullerene& G,
    std::size_t up_to,
    int max_size_l,
    int max_param_sum_b,
    int min_reduction_size)
{
    if (max_size_l < 0) {
        return;
    }

    if (G.total_nodes() >= up_to) {
        return;
    }

    std::vector<std::unique_ptr<base_expansion>> expansions;
    expansions.reserve(512);

    for (int s = 0; s <= max_size_l; s++) {
        {
            auto v = find_l_expansions(G, s);
            expansions.reserve(expansions.size() + v.size());
            for (auto& up : v) {
                expansions.push_back(std::move(up));
            }
        }
    }
    for (int s=0;s<=max_param_sum_b;s++)
    {
            for (int pre = 0; pre <= s; pre++) {
                int post = s - pre;
                auto v = find_b_expansions(G, pre, post);
                expansions.reserve(expansions.size() + v.size());
                for (auto& up : v) {
                    expansions.push_back(std::move(up));
                }
            }
    }



    for (auto& up : expansions) {
        if (!up->validate()) {
            continue;
        }


        if (auto* le = dynamic_cast<l_expansion*>(up.get())) {


            const auto cand = le->candidate();
            up->apply();;
            auto red = std::make_unique<l_reduction>(matching_reduction_from_expansion(*le));

            if (red->is_canonical(G, min_reduction_size)) {
                register_and_emit(G);
                int next_max_l_bound = bound_by_vertex_count_l(G, up_to);
                int next_max_b_bound = bound_by_vertex_count_b(G, up_to);
                int next_max_l = std::min(next_max_l_bound, 4);
                int next_max_b = std::min(next_max_b_bound, 4);
                dfs_(G, up_to, next_max_l, next_max_b, min_reduction_size);
                G.reduce_id();
            }


            red->apply(G, cand);
            continue;
        }

        if (auto* be = dynamic_cast<b_expansion*>(up.get())) {
            const auto cand = be->candidate();
            up->apply();
     
            auto red = std::make_unique<b_reduction>(matching_reduction_from_expansion(*be));
            if (red->is_canonical(G, min_reduction_size)) {
                register_and_emit(G);
                int next_max_l_bound = bound_by_vertex_count_l(G, up_to);
                int next_max_b_bound = bound_by_vertex_count_b(G, up_to);
                int next_max_l = std::min(next_max_l_bound, 4);
                int next_max_b = std::min(next_max_b_bound, 4);
                dfs_(G, up_to, next_max_l, next_max_b, min_reduction_size);
                G.reduce_id();
            }
            
            red->apply(G, cand);
            continue;
        }
    }
    
}
