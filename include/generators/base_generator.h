#ifndef BASE_GENERATOR_H
#define BASE_GENERATOR_H
#include <cstddef>
#include <fullerene/dual_fullerene.h>
#include <iostream>
#include <map>

class base_generator {
public:
    std::map<size_t, int> counts;
    virtual ~base_generator() = default;
    virtual void generate(std::size_t up_to) = 0;
    virtual void register_and_emit(dual_fullerene& G) {
        G.register_id();
        auto P = G.to_primal();
        counts[P.get_size()]++;
        //std::cout << P << std::flush;
    }
};

#endif //BASE_GENERATOR_H
