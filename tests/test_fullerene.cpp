#include <map>
#include <catch2/catch_test_macros.hpp>
#include <catch2/internal/catch_preprocessor_internal_stringify.hpp>
#include <catch2/internal/catch_assertion_handler.hpp>
#include <catch2/internal/catch_decomposer.hpp>
#include <fullerene/construct.h>
#include <fullerene/dual_fullerene.h>

void validate_dual_fullerene(const dual_fullerene& f) {
    REQUIRE(f.get_nodes_5().size() == 12);

    for (auto& u: f.get_nodes_5()) {
        REQUIRE(u->degree() == 5);

        for (auto& weak_v : u->neighbors()) {
            auto v = weak_v.lock();
            REQUIRE(v);

            size_t count = 0;
            for (auto& weak_v2 : u->neighbors())
                if (weak_v2.lock() == v) count++;

            REQUIRE(count == 1);
        }
    }

    for (auto& u: f.get_nodes_6()) {
        REQUIRE(u->degree() == 6);

        for (auto& weak_v : u->neighbors()) {
            auto v = weak_v.lock();
            REQUIRE(v);

            size_t count = 0;
            for (auto& weak_v2 : u->neighbors())
                if (weak_v2.lock() == v) count++;

            REQUIRE(count == 1);
        }
    }
}

void validate_fullerene(const fullerene& f, const dual_fullerene& d) {
    const auto& faces = f.get_adjacency();
    const std::size_t V = d.total_nodes();
    const std::size_t E = (5 * 12 + 6 * (V - 12)) / 2;
    const std::size_t F = E - V + 2;

    REQUIRE(faces.size() == F);

    auto normalize = [](unsigned a, unsigned b) {
        return std::minmax(a, b);
    };

    std::map<std::pair<unsigned, unsigned>, unsigned> edge_count;

    for (unsigned u = 0; u < faces.size(); ++u) {
        const auto& adj = faces[u];

        for (unsigned v : adj) {
            REQUIRE(std::find(faces[v].begin(), faces[v].end(), u) != faces[v].end());
        }

        const unsigned a = adj[0];
        const unsigned b = adj[1];
        const unsigned c = adj[2];

        REQUIRE(a != u);
        REQUIRE(b != u);
        REQUIRE(c != u);

        REQUIRE(a < F);
        REQUIRE(b < F);
        REQUIRE(c < F);

        REQUIRE(a != b);
        REQUIRE(b != c);
        REQUIRE(c != a);

        auto e1 = normalize(u, a);
        auto e2 = normalize(u, b);
        auto e3 = normalize(u, c);

        edge_count[e1]++;
        edge_count[e2]++;
        edge_count[e3]++;
    }
}


// construct tests
TEST_CASE("Base dual fullerenes are structurally valid") {
    const auto d1 = create_c20_fullerene();
    const auto d2 = create_c28_fullerene();
    const auto d3 = create_c30_fullerene();

    validate_dual_fullerene(d1);
    validate_dual_fullerene(d2);
    validate_dual_fullerene(d3);
}

// base_node tests
TEST_CASE("base_node basic functionality") {
    const auto n5 = node_5::create(1);
    const auto n6 = node_6::create(2);

    REQUIRE(n5->expected_degree() == 5);
    REQUIRE(n6->expected_degree() == 6);

    REQUIRE(n5->degree() == 0);
    REQUIRE(n6->degree() == 0);
}

TEST_CASE("Neighbors can be added and are reciprocal") {
    const auto a = node_5::create(1);
    const auto b = node_6::create(2);

    a->add_neighbor(b);
    b->add_neighbor(a);

    REQUIRE(a->degree() == 1);
    REQUIRE(b->degree() == 1);
    REQUIRE(a->is_neighbor_of(b));
    REQUIRE(b->is_neighbor_of(a));
}

// directed_edge tests
TEST_CASE("directed_edge basic navigation") {
    auto a = node_5::create(1);
    auto b = node_6::create(2);

    a->add_neighbor(b);
    b->add_neighbor(a);

    const auto e = a->get_edge(0);

    REQUIRE(e.to() == b);

    const auto inv = e.inverse();
    REQUIRE(inv.to() == a);
}

TEST_CASE("directed_edge complex navigation") {
    auto u = node_5::create_sized(1);
    auto u_neighbors = std::vector<std::shared_ptr<node_5>>(5);

    for (int i = 0; i < u_neighbors.size(); i++) {
        u_neighbors[i] = node_5::create(i);

        u->set_neighbor_at(i, u_neighbors[i]);
        u_neighbors[i]->add_neighbor(u);
    }

    for (int i = 0; i < u_neighbors.size(); i++) {
        u_neighbors[i]->add_neighbor(u_neighbors[(i + 1) % u_neighbors.size()]);
        u_neighbors[i]->add_neighbor(u_neighbors[(i - 1 + u_neighbors.size()) % u_neighbors.size()]);
    }

    auto e = u->get_edge(0);

    REQUIRE(e.to() == u_neighbors[0]);
    REQUIRE(e.inverse().to() == u);
    REQUIRE(e.left_turn(1).to() == u_neighbors[1]);
    REQUIRE(e.right_turn(1).to() == u_neighbors[4]);
    REQUIRE(e.left_turn(1).left_turn(1).to() == u);
    REQUIRE(e.right_turn(1).right_turn(1).to() == u);
    REQUIRE(e.prev_around(1).to() == u_neighbors[4]);
    REQUIRE(e.next_around(1).to() == u_neighbors[1]);

    for (int i = 0; i < u_neighbors.size(); i++) {
        REQUIRE(e.next_around(i).prev_around(i).to() == u_neighbors[0]);
        REQUIRE(e.prev_around(i).next_around(i).to() == u_neighbors[0]);
    }
}

// dual_fullerene tests

// fullerene tests
TEST_CASE("Fullerenes generated from base dual fullerenes are structurally valid") {
    const auto d1 = create_c20_fullerene();
    const auto d2 = create_c28_fullerene();
    const auto d3 = create_c30_fullerene();

    const auto f1 = d1.to_primal();
    const auto f2 = d2.to_primal();
    const auto f3 = d3.to_primal();

    validate_fullerene(f1, d1);
    validate_fullerene(f2, d2);
    validate_fullerene(f3, d3);
}
