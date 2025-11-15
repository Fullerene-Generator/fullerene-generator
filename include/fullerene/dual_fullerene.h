#ifndef DUAL_FULLERENE_H
#define DUAL_FULLERENE_H
#include <fullerene/base_node.h>
#include <memory>
#include <vector>

class dual_fullerene {
    std::vector<std::shared_ptr<node_5>> nodes_5;
    std::vector<std::shared_ptr<node_6>> nodes_6;

public:
    explicit dual_fullerene(const std::vector<std::vector<unsigned int>>& adjacency);

    [[nodiscard]] const std::vector<std::shared_ptr<node_5>>& get_nodes_5() const noexcept { return nodes_5; }
    [[nodiscard]] const std::vector<std::shared_ptr<node_6>>& get_nodes_6() const noexcept { return nodes_6; }
    [[nodiscard]] std::size_t total_nodes() const noexcept { return nodes_5.size() + nodes_6.size(); }
};

#endif //DUAL_FULLERENE_H
