#ifndef MAIN_GENERATOR_H
#define MAIN_GENERATOR_H

#include <generators/base_generator.h>
#include <fullerene/dual_fullerene.h>

class main_generator final : base_generator {
public:
    main_generator() = default;

    void generate(std::size_t up_to) override;

private:
    void dfs_(dual_fullerene& G,
        std::size_t up_to,
        int max_param_sum,
        int min_reduction_size);

    static int bound_by_vertex_count_(const dual_fullerene& G, std::size_t up_to);
};

#endif // MAIN_GENERATOR_H
