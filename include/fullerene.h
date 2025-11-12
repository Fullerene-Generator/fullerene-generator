#ifndef FULLERENE_H
#define FULLERENE_H
#include <memory>
#include <vector>

enum class node_type {
    NODE_5,
    NODE_6,
  };

class base_node {
protected:
    std::vector<std::weak_ptr<base_node>> neighbors_;
    unsigned int id_ = 0;

    explicit base_node(const unsigned int id) : id_(id) {}
public:
    virtual ~base_node() = default;

    [[nodiscard]] virtual std::size_t expected_degree() const = 0;

    [[nodiscard]] unsigned int id() const { return id_; }

    void add_neighbor(const std::shared_ptr<base_node>& n);

    [[nodiscard]] std::size_t neighbor_count() const;
};

class node_5 final : public base_node {
    explicit node_5(const unsigned int id) : base_node(id) {}
public:
    static std::shared_ptr<node_5> create(const unsigned int id) {
        return std::shared_ptr<node_5>(new node_5(id));
    }

    [[nodiscard]] std::size_t expected_degree() const override { return 5; }
};

class node_6 final : public base_node {
    explicit node_6(const unsigned int id) : base_node(id) {}
public:
    static std::shared_ptr<node_6> create(const unsigned int id) {
        return std::shared_ptr<node_6>(new node_6(id));
    }

    [[nodiscard]] std::size_t expected_degree() const override { return 6; }
};

class dual_fullerene {
    std::vector<std::shared_ptr<node_5>> nodes_5;
    std::vector<std::shared_ptr<node_6>> nodes_6;

public:
    explicit dual_fullerene(const std::vector<std::vector<unsigned int>>& adjacency);

    [[nodiscard]] const std::vector<std::shared_ptr<node_5>>& get_nodes_5() const noexcept;
    [[nodiscard]] const std::vector<std::shared_ptr<node_6>>& get_nodes_6() const noexcept;

    [[nodiscard]] std::size_t total_nodes() const noexcept { return nodes_5.size() + nodes_6.size(); }
};

#endif //FULLERENE_H
