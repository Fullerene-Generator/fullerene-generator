#ifndef CONSTRUCT_H
#define CONSTRUCT_H

#include <fullerene.h>
#include <vector>

const std::vector<std::vector<unsigned int>> dodecahedron_adjacency = {
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
    return dual_fullerene(dodecahedron_adjacency);
}

#endif //CONSTRUCT_H
