#ifndef FULLERENE_GENERATOR_ID_REGISTRY_H
#define FULLERENE_GENERATOR_ID_REGISTRY_H
#include <map>
#include <string>
#include <sstream>

class id_registry {
    inline static std::map<unsigned, unsigned> counters;

public:
    [[nodiscard]] static std::string register_id(unsigned n) {
        const auto id = counters[n]++;

        std::ostringstream ss;

        ss << n << ":" << id;
        return ss.str();
    }
};

#endif //FULLERENE_GENERATOR_ID_REGISTRY_H
