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

} // namespace

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
        dfs_(G, up_to, 4, 4, 1);
    }

    for (const auto& [v, c] : counts_) {
        std::cout << v << " vertices: " << c << "\n";
    }
}

int main_generator::bound_by_vertex_count_l(const dual_fullerene& G, std::size_t up_to)
{
    // primal vertices = 20 + 2 * (#hexagons in dual)
    int primal_v = 20 + 2 * static_cast<int>(G.get_nodes_6().size());
    int dif = static_cast<int>(up_to) - primal_v;

    // L_i adds 2*(i+2) primal vertices, B_{p,q} adds 2*(p+q+2)
    return dif / 2 - 2;
}

int main_generator::bound_by_vertex_count_b(const dual_fullerene& G, std::size_t up_to)
{
    // primal vertices = 20 + 2 * (#hexagons in dual)
    int primal_v = 20 + 2 * static_cast<int>(G.get_nodes_6().size());
    int dif = static_cast<int>(up_to) - primal_v;

    // L_i adds 2*(i+2) primal vertices, B_{p,q} adds 2*(p+q+2)
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

    //std::cout << "found " << expansions.size() << " expansions\n";

    int l_count = 0;
    int l_can_count = 0;
    int b_count = 0;
    int b_can_count = 0;

    for (auto& up : expansions) {
        if (!up->validate()) {
            continue;
        }
        //std::cout << "validated\n";

        // ---- L expansion ----
        if (auto* le = dynamic_cast<l_expansion*>(up.get())) {
           // std::cout << "l expansion size: " << le->candidate().length << '\n';
          //  std::cout << "expansion first edge, from: " << le->candidate().start.from->id()
            //    << " to: " << le->candidate().start.to()->id() << '\n';

            // 1) copy candidate BEFORE apply (safe even if apply mutates internal state)
            l_count++;
            const auto cand = le->candidate();

            // 2) apply expansion (this is what initializes inverse_*_edge)
            up->apply();
            //std::cout << "expansion applied\n";

            // 3) now it is safe to read inverse edges + construct reduction
            auto inv1 = le->inverse_first_edge();
            auto inv2 = le->inverse_second_edge();
            //std::cout << "expansion inverse first edge from: " << inv1.from->id()
              //  << " to: " << inv1.to()->id()
                //<< " second edge " << inv2.from->id()
                //<< " to " << inv2.to()->id() << '\n';

            auto red = std::make_unique<l_reduction>(matching_reduction_from_expansion(*le));

            //std::cout << "min reduction size: " << min_reduction_size << '\n';
            if (red->is_canonical(G, min_reduction_size)) {
              //  std::cout << "red is canonical\n";
                l_can_count++;
                register_and_emit(G);
                int next_max_l_bound = bound_by_vertex_count_l(G, up_to);
                int next_max_b_bound = bound_by_vertex_count_b(G, up_to);
                int next_max_l = std::min(next_max_l_bound, 8);
                int next_max_b = std::min(next_max_l_bound, 8);
                dfs_(G, up_to, next_max_l, next_max_b, min_reduction_size);
                G.reduce_id();
            }


           // std::cout << "tries to apply\n";
            red->apply(G, cand);
            //std::cout << "reduction applied\n";
            continue;
        }

        // ---- B expansion ----
        if (auto* be = dynamic_cast<b_expansion*>(up.get())) {
           // std::cout << "b expansion size: " << be->candidate().length_pre_bend
             //   << ' ' << be->candidate().length_post_bend << '\n';

            // 1) copy candidate BEFORE apply
            b_count++;
            const auto cand = be->candidate();

            // 2) apply expansion (initializes inverse edges for B too)
            up->apply();
            //std::cout << "expansion applied\n";

            // 3) construct reduction AFTER apply
            auto red = std::make_unique<b_reduction>(matching_reduction_from_expansion(*be));

            //std::cout << "min reduction size: " << min_reduction_size << '\n';
            if (red->is_canonical(G, min_reduction_size)) {
              //  std::cout << "red is canonical\n";
                b_can_count++;
                register_and_emit(G);
                int next_max_l_bound = bound_by_vertex_count_l(G, up_to);
                int next_max_b_bound = bound_by_vertex_count_b(G, up_to);
                int next_max_l = std::min(next_max_l_bound, 4);
                int next_max_b = std::min(next_max_l_bound, 8);
                dfs_(G, up_to, next_max_l, next_max_b, min_reduction_size);
                G.reduce_id();
            }
            

           // std::cout << "tries to apply\n";
            red->apply(G, cand);
            //std::cout << "reduction applied\n";
            continue;
        }
    }
    if (true || b_can_count > 0 || (G.get_nodes_6().size()*2 >= 6 && G.get_nodes_6().size() <= 6)) {
        std::cout << "graph size: " << G.get_nodes_6().size() * 2 + 20 << " l count: " << l_count << " l can count: " << l_can_count << " b count: " << b_count << " b can count: " << b_can_count << '\n';
    }
    
}
