#include <fullerene/base_node.h>
#include <fullerene/directed_edge.h>

std::shared_ptr<base_node> directed_edge::to() const {
    return from->neighbor_at(index);
}

directed_edge directed_edge::inverse() const {
    return  to()->get_edge(from);
}

directed_edge directed_edge::next_around(const unsigned int times) const {
    return { from, (index + times) % from->degree() };
}

directed_edge directed_edge::prev_around(const unsigned int times) const {
    const std::size_t d = from->degree();
    return { from, (index - times + d) % d };
}

directed_edge directed_edge::left_turn(const unsigned int which) const {
    return inverse().next_around(which);
}

directed_edge directed_edge::right_turn(const unsigned int which) const {
    return inverse().prev_around(which);
}

edge_data &directed_edge::data() const {
    return from->get_edge_data(index);
}
