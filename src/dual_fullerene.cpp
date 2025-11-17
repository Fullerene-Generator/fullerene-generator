#include <fullerene/dual_fullerene.h>
#include <string>
#include <stdexcept>
#include <algorithm>

template<typename F>
void dual_fullerene::for_each_node(F&& f) const {
    for (const auto& node : nodes_5) f(node);
    for (const auto& node : nodes_6) f(node);
}

dual_fullerene::dual_fullerene(const std::vector<std::vector<unsigned int>>& adjacency) {
    const std::size_t n = adjacency.size();

    std::vector<std::shared_ptr<base_node>> index_to_node(n);

    for (std::size_t i = 0; i < n; i++) {
        if (i < 12) {
            auto p = node_5::create(static_cast<unsigned int>(i));
            index_to_node[i] = p;
            nodes_5.push_back(p);
        }
        else {
            auto p = node_6::create(static_cast<unsigned int>(i));
            index_to_node[i] = p;
            nodes_6.push_back(p);
        }
    }

    for (std::size_t i = 0; i < n; i++) {
        const auto& neighs = adjacency[i];

        for (const auto j : neighs) {
            if (j >= n)
                throw std::out_of_range("Adjacency index " + std::to_string(j) + " out of range");
            if (j == i)
                throw std::invalid_argument("Self-loop at node " + std::to_string(i));
        }

        const std::size_t deg = neighs.size();
        if (i < 12 && deg != 5)
            throw std::invalid_argument("5-node " + std::to_string(i) + " has degree " +
                std::to_string(deg) + ", expected 5");
        if (i >= 12 && deg != 6)
            throw std::invalid_argument("6-node " + std::to_string(i) + " has degree " +
                std::to_string(deg) + ", expected 6");
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (auto j : adjacency[i]) {
            const auto& neighs = adjacency[j];
            if (std::ranges::find(neighs.begin(), neighs.end(), static_cast<unsigned int>(i)) == neighs.end()) {
                throw std::invalid_argument("Adjacency not symmetric between " +
                    std::to_string(i) + " and " + std::to_string(j));
            }
        }
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (const auto j : adjacency[i]) {
            const auto& a = index_to_node[i];
            const auto& b = index_to_node[j];
            a->add_neighbor(b);
        }
    }
}

fullerene dual_fullerene::to_primal() const {
    const std::size_t V = total_nodes();
    const std::size_t E = (5 * 12 + 6 * (V - 12)) / 2;
    const std::size_t F = E - V + 2;

    std::vector<std::array<unsigned int, 3>> adjacency(F);
    std::vector<unsigned int> counts(F, 0);
    unsigned int face = 0;

    for_each_node([&](const std::shared_ptr<base_node>& node) {
        for (int i = 0; i < static_cast<int>(node->degree()); i++) {
            auto edge = node->get_edge(static_cast<std::size_t>(i));
            if (edge.data().marked) continue;

            const auto start_node = edge.from;

            do {
                edge.data().marked = true;
                edge.data().rhs_face_index = face;
                edge = edge.right_turn();
            } while (edge.from != start_node);

            face++;
        }
        });

    for_each_node([&](const std::shared_ptr<base_node>& node) {
        for (int i = 0; i < static_cast<int>(node->degree()); i++) {
            auto edge = node->get_edge(static_cast<std::size_t>(i));

            const auto u = edge.data().rhs_face_index;
            const auto v = edge.inverse().data().rhs_face_index;

            adjacency[u][counts[u]++] = v;
        }
        });

    const auto outer_face = nodes_5[0];
    auto outer_face_nodes = std::array<unsigned int, 5>();

    for (int i = 0; i < static_cast<int>(outer_face->degree()); i++) {
        auto edge = outer_face->get_edge(static_cast<std::size_t>(i));
        outer_face_nodes[i] = edge.data().rhs_face_index;
    }

    clear_all_edge_data();

    return fullerene(adjacency, outer_face_nodes);
}

std::shared_ptr<base_node> dual_fullerene::get_node(unsigned int id) const {
    for (const auto& n : nodes_5) {
        if (n->id() == id)
            return n;
    }
    for (const auto& n : nodes_6) {
        if (n->id() == id)
            return n;
    }

    throw std::out_of_range("No node with id " + std::to_string(id));
}

void dual_fullerene::clear_all_edge_data() const {
    for_each_node([](const std::shared_ptr<base_node>& node) {
        node->clear_all_edge_data();
        });
}

int dual_fullerene::add_vertex(node_type type) {
    const unsigned int id = static_cast<unsigned int>(total_nodes());

    switch (type) {
    case node_type::NODE_5: {
        auto p = node_5::create(id);
        nodes_5.push_back(p);
        break;
    }
    case node_type::NODE_6: {
        auto p = node_6::create(id);
        nodes_6.push_back(p);
        break;
    }
    }

    return static_cast<int>(id);
}


void dual_fullerene::add_neighbour_after(int v, int after, int v2) {
    auto v_node = get_node(static_cast<unsigned int>(v));
    auto after_node = get_node(static_cast<unsigned int>(after));
    auto v2_node = get_node(static_cast<unsigned int>(v2));
    v_node->add_neighbour_after(after_node, v2_node);
}

void dual_fullerene::add_neighbour_before(int v, int before, int v2) {
    auto v_node = get_node(static_cast<unsigned int>(v));
    auto before_node = get_node(static_cast<unsigned int>(before));
    auto v2_node = get_node(static_cast<unsigned int>(v2));
    v_node->add_neighbour_before(before_node, v2_node);
}

void dual_fullerene::remove_edge(int v1, int v2) {
    auto n1 = get_node(static_cast<unsigned int>(v1));
    auto n2 = get_node(static_cast<unsigned int>(v2));
    n1->remove_neighbor(n2);
    n2->remove_neighbor(n1);
}

void dual_fullerene::replace_neighbour(int v, int old_n, int new_n) {
    auto v_node = get_node(static_cast<unsigned int>(v));
    auto old_node = get_node(static_cast<unsigned int>(old_n));
    auto new_node = get_node(static_cast<unsigned int>(new_n));
    v_node->replace_neighbor(old_node, new_node);
}

void dual_fullerene::move_neighbourhood(int from, int to) {
    auto from_node = get_node(static_cast<unsigned int>(from));
    auto to_node = get_node(static_cast<unsigned int>(to));

    to_node->move_neighborhood_from(from_node);

    for (const auto& w : to_node->neighbors()) {
        if (auto n = w.lock()) {
            n->replace_neighbor(from_node, to_node);
        }
    }
}
