#include "fullerene.h"

#include <utility>

dual_fullerene::dual_fullerene(std::vector<std::shared_ptr<node>> nodes) {
    this->nodes = std::move(nodes);
}

void dual_fullerene::replace_neighbor(int u, int oldN, int newN)
{
}

void dual_fullerene::insert_neighbor_after(int u, int afterN, int newN)
{
}

void dual_fullerene::add_edge_between(int x, int x_after, int y, int y_after)
{
}

int dual_fullerene::add_vertex(node_type t)
{
    return 0;
}
