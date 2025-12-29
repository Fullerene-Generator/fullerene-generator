#ifndef BASE_GENERATOR_H
#define BASE_GENERATOR_H
#include <cstddef>
#include <fullerene/dual_fullerene.h>
#include <iostream>

class base_generator {
public:
    virtual ~base_generator() = default;
    virtual void generate(std::size_t up_to) = 0;
    [[nodiscard]] static std::string register_and_emit(const dual_fullerene& G, const std::string& parent_id = "BASE") {
        auto P = G.to_primal();
        auto id = P.register_fullerene(parent_id);
        std::cout << P << std::flush;

        return id;
    }
};

#endif //BASE_GENERATOR_H
