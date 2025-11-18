#ifndef BASE_NODE_H
#define BASE_NODE_H
#include <fullerene/directed_edge.h>
#include <memory>
#include <vector>

enum class node_type {
    NODE_5,
    NODE_6,
  };

class base_node : public std::enable_shared_from_this<base_node> {
protected:
    std::vector<std::weak_ptr<base_node>> neighbors_;
    std::vector<edge_data> edges_;
    unsigned int id_ = 0;
    node_type type_;

    explicit base_node(const unsigned int id, const node_type type, const std::size_t neighbors) : id_(id), type_(type) { resize_(neighbors); }

    void resize_(std::size_t neighbors);

public:
    virtual ~base_node() = default;
    [[nodiscard]] virtual std::size_t expected_degree() const = 0;
    [[nodiscard]] unsigned int id() const { return id_; }
    [[nodiscard]] node_type type() const { return type_; }
    [[nodiscard]] std::size_t degree() const;
    [[nodiscard]] std::shared_ptr<base_node> neighbor_at(std::size_t index) const;
    [[nodiscard]] const std::vector<std::weak_ptr<base_node>>& neighbors() const;
    [[nodiscard]] directed_edge get_edge(std::size_t index);
    [[nodiscard]] directed_edge get_edge(const std::shared_ptr<const base_node>& other);
    [[nodiscard]] bool is_neighbor_of(const std::shared_ptr<base_node> &other) const;
    void set_neighbor_at(std::size_t index, const std::shared_ptr<base_node>& n);
    [[nodiscard]] edge_data& get_edge_data(std::size_t index);
    void clear_all_edge_data();
};

class node_5 final : public base_node {
    explicit node_5(const unsigned int id) : base_node(id, node_type::NODE_5, 5) {}
public:
    static std::shared_ptr<node_5> create(const unsigned int id) {
        return std::shared_ptr<node_5>(new node_5(id));
    }

    [[nodiscard]] std::size_t expected_degree() const override { return 5; }
};

class node_6 final : public base_node {
    explicit node_6(const unsigned int id) : base_node(id, node_type::NODE_6, 6) {}
public:
    static std::shared_ptr<node_6> create(const unsigned int id) {
        return std::shared_ptr<node_6>(new node_6(id));
    }

    [[nodiscard]] std::size_t expected_degree() const override { return 6; }
};

#endif //BASE_NODE_H
