#ifndef FULLERENE_GENERATOR_B_EXPANSION_H
#define FULLERENE_GENERATOR_B_EXPANSION_H
#include "base_expansion.h"


class b_expansion: base_expansion {
public:
    explicit b_expansion(dual_fullerene &G)
        : base_expansion(G) {
    }

    [[nodiscard]] bool validate() const override;
    void apply() override;
};


#endif //FULLERENE_GENERATOR_B_EXPANSION_H