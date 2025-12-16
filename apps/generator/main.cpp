#include "generators/f_expansion_generator.h"
#include <iostream>
#include <string>
#include "generators//l_expansion_generator.h"


int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: fullerene_generator <max_size>\n";
        return 1;
    }

    size_t max_size = std::stoul(argv[1]);

    auto generator = f_expansion_generator();
    generator.generate(max_size);

    auto generator_l = l_expansion_generator();
    generator_l.generate(max_size);
}
