#include <ranges>
#include <sstream>
#include <fullerene/fullerene.h>

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
