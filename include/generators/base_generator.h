#ifndef BASE_GENERATOR_H
#define BASE_GENERATOR_H
#include <cstddef>

class base_generator {
public:
    virtual ~base_generator() = default;
    virtual void generate(std::size_t up_to) = 0;
};

#endif //BASE_GENERATOR_H
