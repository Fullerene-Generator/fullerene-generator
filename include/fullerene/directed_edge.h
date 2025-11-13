#ifndef DIRECTED_EDGE_H
#define DIRECTED_EDGE_H
#include <memory>

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

#endif //DIRECTED_EDGE_H
