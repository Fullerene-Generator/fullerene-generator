#ifndef BASE_EXPANSION_H
#define BASE_EXPANSION_H
#include <fullerene/dual_fullerene.h>

bool patch_nodes_unique(const std::vector<int>& path, const std::vector<int>& para);

class base_expansion {
protected:
    dual_fullerene& G_;

    explicit base_expansion(dual_fullerene& G) : G_(G) {}

public:
    virtual ~base_expansion() = default;

    [[nodiscard]] virtual bool validate() const = 0;
    virtual void apply() = 0;
};

#endif //BASE_EXPANSION_H
