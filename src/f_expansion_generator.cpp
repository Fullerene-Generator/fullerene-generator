#include <iostream>
#include <expansions/f_expansion.h>
#include <fullerene/construct.h>
#include <generator/f_expansion_generator.h>

constexpr size_t C30_SIZE = 20;
constexpr size_t F_EXPANSION_SIZE_INCREMENT = 10;

void f_expansion_generator::generate(std::size_t up_to) {
    if (up_to < C30_SIZE) {
        return;
    }

    size_t current_size = C30_SIZE;
    auto C30_dual = create_c30_fullerene();

    auto C30 = C30_dual.to_primal();
    C30.compute_tutte_embedding();
    std::cout << C30 << std::flush;

    while (current_size < up_to) {
        auto expansion = f_expansion(C30_dual, C30_dual.get_nodes_5()[0]);
        if (!expansion.validate()) {
            throw std::logic_error("The F expansion can't be performed");
        }
        expansion.apply();

        C30 = C30_dual.to_primal();
        C30.compute_tutte_embedding();
        std::cout << C30;

        current_size += F_EXPANSION_SIZE_INCREMENT;
    }
}
