#include <fullerene.h>
#include <vector>

constexpr unsigned int C20_dual_fullerene_vertex_count = 12;

const unsigned int dodecahedron_adjacency[12][5] = {
    {1, 4, 5, 2, 3},
    {0, 3, 7, 6, 4},
    {0, 5, 9, 8, 3},
    {0, 2, 8, 7, 1},
    {0, 1, 6, 10, 5},
    {0, 4, 10, 9, 2},
    {1, 7, 11, 10, 4},
    {1, 3, 8, 11, 6},
    {2, 9, 11, 7, 3},
    {2, 5, 10, 11, 8},
    {4, 6, 11, 9, 5},
    {6, 7, 8, 9, 10}
};

dual_fullerene create_C20_fullerene() {
    auto nodes = std::vector<std::shared_ptr<node>>(C20_dual_fullerene_vertex_count);

    for (unsigned int i = 0; i < C20_dual_fullerene_vertex_count; i++) {
        nodes[i] = std::make_shared<node>(node {
            .type = node_type::NODE_5,
            .neighbors = std::vector<std::shared_ptr<node>>(5),
        });
    }

    for (unsigned int i = 0; i < C20_dual_fullerene_vertex_count; i++) {
        for (unsigned int j = 0; j < 5; j++) {
            const unsigned int neighbor_index = dodecahedron_adjacency[i][j];
            nodes[i]->neighbors[j] = nodes[neighbor_index];
        }
    }

    return dual_fullerene(nodes);
}
