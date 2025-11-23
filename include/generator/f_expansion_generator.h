#ifndef F_EXPANSION_GENERATOR_H
#define F_EXPANSION_GENERATOR_H
#include <generator/base_generator.h>

class f_expansion_generator final : base_generator {
public:
    f_expansion_generator() = default;

    void generate(std::size_t up_to) override;
};

#endif //F_EXPANSION_GENERATOR_H
