#ifndef L_EXPANSION_GENERATOR_H
#define L_EXPANSION_GENERATOR_H

#include <iostream>
#include <generators/base_generator.h>
#include <fullerene/dual_fullerene.h>

class l_expansion_generator final : base_generator {
public:
    l_expansion_generator() = default;

    void generate(std::size_t up_to) override;

private:
    void dfs_(dual_fullerene& G,
        std::size_t up_to,
        int max_i,
        int min_reduction_size);

    static int bound_by_vertex_count_(const dual_fullerene& G, std::size_t up_to);
};

#endif // L_EXPANSION_GENERATOR_H
