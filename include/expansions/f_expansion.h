#ifndef F_EXPANSION_H
#define F_EXPANSION_H
#include <memory>
#include <utility>
#include <expansions/base_expansion.h>
#include <fullerene/base_node.h>

class f_expansion final : public base_expansion {
    std::shared_ptr<node_5> v_;

public:
    explicit f_expansion(dual_fullerene &G,  std::shared_ptr<node_5> v): base_expansion(G), v_(std::move(v)) {}

    [[nodiscard]] bool validate() const override;
    void apply() override;
};

#endif //F_EXPANSION_H
