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
        int max_param_sum = bound_by_vertex_count_(G, up_to);
        int max0 = std::min(max_param_sum, 3);
        dfs_(G, up_to, max0, 1);
    }

    if (up_to < 28) {
        return;
    }

    {
        auto G = create_c28_fullerene();
        register_and_emit(G);
        int max_param_sum = bound_by_vertex_count_(G, up_to);
        int max0 = std::min(max_param_sum, 4);
        dfs_(G, up_to, max0, 1);
    }
}

int main_generator::bound_by_vertex_count_(const dual_fullerene& G, std::size_t up_to)
{
    // primal vertices = 20 + 2 * (#hexagons in dual)
    int primal_v = 20 + 2 * static_cast<int>(G.get_nodes_6().size());
    int dif = static_cast<int>(up_to) - primal_v;

    // L_i adds 2*(i+2) primal vertices, B_{p,q} adds 2*(p+q+2)
    return dif / 2 - 2;
}

void main_generator::dfs_(dual_fullerene& G,
    std::size_t up_to,
    int max_param_sum,
    int min_reduction_size)
{
    if (max_param_sum < 0) {
        return;
    }

    std::vector<std::unique_ptr<base_expansion>> expansions;
    expansions.reserve(512);

    // For a fixed size s:
    //   - all L expansions with length = s
    //   - all B expansions with pre+post = s
    for (int s = 0; s <= max_param_sum; s++) {
        {
            auto v = find_l_expansions(G, s);
            expansions.reserve(expansions.size() + v.size());
            for (auto& up : v) {
                expansions.push_back(std::move(up));
            }
        }
        if (s == 0) continue;

        {
            for (int pre = 0; pre <= s-1; pre++) {
                int post = s - 1 - pre;
                auto v = find_b_expansions(G, pre, post);
                expansions.reserve(expansions.size() + v.size());
                for (auto& up : v) {
                    expansions.push_back(std::move(up));
                }
            }
        }
    }

    //std::cout << "found " << expansions.size() << " expansions\n";

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
                register_and_emit(G);
                int next_max = bound_by_vertex_count_(G, up_to);
                int next_max0 = std::min(next_max, 8);
                dfs_(G, up_to, next_max0, min_reduction_size);
                G.reduce_id();
            }


           // std::cout << "tries to apply\n";
            red->apply(G, cand);
            //std::cout << "reduction applied\n";
            continue;
        }

        // ---- B expansion ----
        if (auto* be = dynamic_cast<b_expansion*>(up.get())) {
            std::cout << "b expansion size: " << be->candidate().length_pre_bend
                << ' ' << be->candidate().length_post_bend << '\n';

            // 1) copy candidate BEFORE apply
            const auto cand = be->candidate();

            // 2) apply expansion (initializes inverse edges for B too)
            up->apply();
            //std::cout << "expansion applied\n";

            // 3) construct reduction AFTER apply
            auto red = std::make_unique<b_reduction>(matching_reduction_from_expansion(*be));

            //std::cout << "min reduction size: " << min_reduction_size << '\n';
            if (red->is_canonical(G, min_reduction_size)) {
              //  std::cout << "red is canonical\n";
                register_and_emit(G);
                int next_max = bound_by_vertex_count_(G, up_to);
                int next_max0 = std::min(next_max, 8);
                dfs_(G, up_to, next_max0, min_reduction_size);
                G.reduce_id();
            }
            

           // std::cout << "tries to apply\n";
            red->apply(G, cand);
            //std::cout << "reduction applied\n";
            continue;
        }
    }
}
