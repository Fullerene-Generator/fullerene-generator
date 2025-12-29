#include <ranges>
#include <sstream>
#include <fullerene/fullerene.h>
#include <generators/id_registry.h>

std::string fullerene::write_all() const noexcept {
    std::ostringstream ss;

    // fullerene metadata
    ss << get_size() << " " << id_ << " " << parent_id_ << "\n";

    // nodes of the outer face
    for (unsigned v : outer_face_nodes_) {
        ss << v << ' ';
    }
    ss << "\n";

    // adjacency of vertices
    for (auto const& adj : adjacency_) {
        ss << adj[0] << " " << adj[1] << " " << adj[2] << "\n";
    }

    return ss.str();
}

std::ostream & operator<<(std::ostream &os, const fullerene &f) {
    os << f.write_all();
    return os;
}

std::string fullerene::register_fullerene(const std::string& parent_id) {
    id_ = id_registry::register_id(get_size());
    parent_id_ = parent_id;

    return id_;
}
