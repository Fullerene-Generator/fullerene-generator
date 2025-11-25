#include <cstdlib>
#include <generator/f_expansion_generator.h>


int main(int argc, char** argv) {
    auto generator = f_expansion_generator();

    if (argc > 1) {
        generator.generate(static_cast<size_t>(strtol(argv[1], nullptr, 10)));
    }

    return 0;
}
