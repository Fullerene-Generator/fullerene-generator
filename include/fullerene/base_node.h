#ifndef BASE_NODE_H
#define BASE_NODE_H
#include <fullerene/directed_edge.h>
#include <memory>
#include <optional>
#include <vector>

enum class node_type {
    NODE_5,
    NODE_6,
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

#endif //BASE_NODE_H
