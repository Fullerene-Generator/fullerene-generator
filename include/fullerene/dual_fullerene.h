#ifndef DUAL_FULLERENE_H
#define DUAL_FULLERENE_H
#include <fullerene/base_node.h>
#include <fullerene/fullerene.h>
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
    [[nodiscard]] fullerene to_primal() const;
    [[nodiscard]] std::shared_ptr<base_node> get_node(unsigned int id) const;
    template<typename F>
    void for_each_node(F&& f) const {
        for (const auto& node : nodes_5) f(node);
        for (const auto& node : nodes_6) f(node);
    }
    void clear_all_edge_data() const;
    int add_vertex(node_type type);
    void add_neighbor_after(int v, int after, int v2);
    void add_neighbor_before(int v, int before, int v2);
    void remove_edge(int v1, int v2);
    void replace_neighbor(int v, int old_n, int new_n);
    void move_neighborhood(int from, int to);
    void add_node(const std::shared_ptr<node_6>& new_node);
    void pop_last_node6();
};

#endif //DUAL_FULLERENE_H
