#include <validation.h>
#include <catch2/internal/catch_preprocessor_internal_stringify.hpp>
#include <catch2/internal/catch_test_macro_impl.hpp>
#include <catch2/internal/catch_test_registry.hpp>

// construct tests
TEST_CASE("Base dual fullerenes are structurally valid", "[dual_fullerene]") {
    const auto d1 = create_c20_fullerene();
    const auto d2 = create_c28_fullerene();
    const auto d3 = create_c30_fullerene();

    validate_dual_fullerene(d1);
    validate_dual_fullerene(d2);
    validate_dual_fullerene(d3);
}

// base_node tests
TEST_CASE("base_node basic functionality", "[base_node]") {
    const auto n5 = node_5::create(1);
    const auto n6 = node_6::create(2);

    REQUIRE(n5->expected_degree() == 5);
    REQUIRE(n6->expected_degree() == 6);

    REQUIRE(n5->degree() == 0);
    REQUIRE(n6->degree() == 0);
}

TEST_CASE("Neighbors can be added and are reciprocal", "[base_node]") {
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
TEST_CASE("directed_edge inverse is involutive", "[directed_edge]") {
    auto a = node_5::create(0);
    auto b = node_6::create(1);

    a->add_neighbor(b);
    b->add_neighbor(a);

    auto e = a->get_edge(0);

    REQUIRE(e.inverse().inverse().from == e.from);
    REQUIRE(e.inverse().inverse().to() == e.to());
}

TEST_CASE("directed_edge next/prev around are inverses", "[directed_edge]") {
    auto u = node_5::create_sized(0);
    std::vector<std::shared_ptr<node_5>> nbrs(5);

    for (int i = 0; i < 5; ++i) {
        nbrs[i] = node_5::create(i + 1);
        u->set_neighbor_at(i, nbrs[i]);
        nbrs[i]->add_neighbor(u);
    }

    auto e = u->get_edge(2);

    for (int k = 0; k < 5; ++k) {
        REQUIRE(e.next_around(k).prev_around(k).to() == e.to());
        REQUIRE(e.prev_around(k).next_around(k).to() == e.to());
    }
}

TEST_CASE("directed_edge basic navigation", "[directed_edge]") {
    auto a = node_5::create(1);
    auto b = node_6::create(2);

    a->add_neighbor(b);
    b->add_neighbor(a);

    const auto e = a->get_edge(0);

    REQUIRE(e.to() == b);

    const auto inv = e.inverse();
    REQUIRE(inv.to() == a);
}

TEST_CASE("directed_edge complex navigation", "[directed_edge]") {
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
TEST_CASE("dual_fullerene edges are bidirectional with matching indices", "[dual_fullerene]") {
    const auto d = create_c30_fullerene();

    d.for_each_node([&](const std::shared_ptr<base_node>& u) {
        for (std::size_t i = 0; i < u->degree(); ++i) {
            auto e = u->get_edge(i);
            auto v = e.to();

            REQUIRE(v->is_neighbor_of(u));

            auto inv = e.inverse();
            REQUIRE(inv.to() == u);
        }
    });
}

TEST_CASE("edge_data is per-directed-edge and survives traversal", "[dual_fullerene]") {
    auto d = create_c20_fullerene();

    d.for_each_node([&](const std::shared_ptr<base_node>& u) {
        for (std::size_t i = 0; i < u->degree(); ++i) {
            auto e = u->get_edge(i);
            e.data().marked = true;
        }
    });

    d.for_each_node([&](const std::shared_ptr<base_node>& u) {
        for (std::size_t i = 0; i < u->degree(); ++i) {
            auto e = u->get_edge(i);
            REQUIRE(e.data().marked);
        }
    });

    d.clear_all_edge_data();

    d.for_each_node([&](const std::shared_ptr<base_node>& u) {
        for (std::size_t i = 0; i < u->degree(); ++i) {
            auto e = u->get_edge(i);
            REQUIRE_FALSE(e.data().marked);
        }
    });
}

// fullerene tests
TEST_CASE("Fullerenes generated from base dual fullerenes are structurally valid", "[fullerene]") {
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
