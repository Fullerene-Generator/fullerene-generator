#ifndef FULLERENE_H
#define FULLERENE_H
#include <memory>
#include <optional>
#include <vector>

enum class node_type {
    NODE_5,
    NODE_6,
  };

class base_node;

struct directed_edge {
    std::shared_ptr<base_node> from;
    std::size_t index;

    directed_edge(std::shared_ptr<base_node> from, const std::size_t index) : from(std::move(from)), index(index) {}

    [[nodiscard]] std::shared_ptr<base_node> to() const;
    [[nodiscard]] directed_edge inverse() const;
    [[nodiscard]] directed_edge next_around(unsigned int times = 1) const;
    [[nodiscard]] directed_edge prev_around(unsigned int times = 1) const;
    [[nodiscard]] directed_edge left_turn(unsigned int which = 1) const;
    [[nodiscard]] directed_edge right_turn(unsigned int which = 1) const;
};

class base_node : public std::enable_shared_from_this<base_node> {
protected:
    std::vector<std::weak_ptr<base_node>> neighbors_;
    unsigned int id_ = 0;
    node_type type_;

    explicit base_node(const unsigned int id, const node_type type) : id_(id), type_(type) {}
public:
    virtual ~base_node() = default;
    [[nodiscard]] virtual std::size_t expected_degree() const = 0;
    [[nodiscard]] unsigned int id() const { return id_; }
    [[nodiscard]] node_type type() const { return type_; }
    [[nodiscard]] std::size_t degree() const;
    [[nodiscard]] std::shared_ptr<base_node> neighbor_at(std::size_t index) const;
    [[nodiscard]] const std::vector<std::weak_ptr<base_node>>& neighbors() const;
    [[nodiscard]] std::optional<directed_edge> get_edge(std::size_t index);
    [[nodiscard]] std::optional<directed_edge> get_edge(const std::shared_ptr<const base_node>& neighbor);
    void add_neighbor(const std::shared_ptr<base_node>& n);
};

class node_5 final : public base_node {
    explicit node_5(const unsigned int id) : base_node(id, node_type::NODE_5) {}
public:
    static std::shared_ptr<node_5> create(const unsigned int id) {
        return std::shared_ptr<node_5>(new node_5(id));
    }

    [[nodiscard]] std::size_t expected_degree() const override { return 5; }
};

class node_6 final : public base_node {
    explicit node_6(const unsigned int id) : base_node(id, node_type::NODE_6) {}
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
