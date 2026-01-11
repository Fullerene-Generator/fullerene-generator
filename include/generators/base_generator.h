#ifndef BASE_GENERATOR_H
#define BASE_GENERATOR_H
#include <cstddef>
#include <fullerene/dual_fullerene.h>
#include <iostream>

class base_generator {
public:
    virtual ~base_generator() = default;
    virtual void generate(std::size_t up_to) = 0;
    virtual void register_and_emit(dual_fullerene& G) {
        G.register_id();
        auto P = G.to_primal();
        std::cout << P << std::flush;
    }
};

#endif //BASE_GENERATOR_H
