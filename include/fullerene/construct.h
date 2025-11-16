#ifndef CONSTRUCT_H
#define CONSTRUCT_H
#include <fullerene/dual_fullerene.h>
#include <vector>

const std::vector<std::vector<unsigned int>> c20_adjacency = {
    {1, 2, 3, 4, 5},
    {0, 5, 6, 7, 2},
    {0, 1, 7, 8, 3},
    {0, 2, 8, 9, 4},
    {0, 3, 9, 10, 5},
    {0, 4, 10, 6, 1},
    {1, 5, 10, 11, 7},
    {1, 6, 11, 8, 2},
    {2, 7, 11, 9, 3},
    {3, 8, 11, 10, 4},
    {4, 9, 11, 6, 5},
    {10, 9, 8, 7, 6}
};

const std::vector<std::vector<unsigned int>>  c28_adjacency = {
    {1, 2, 14, 3, 12},
    {0, 12, 4, 13, 2},
    {0, 1, 13, 5, 14},
    {0, 14, 6, 7, 12},
    {1, 12, 8, 9, 13},
    {2, 13, 10, 11, 14},
    {3, 14, 11, 15, 7},
    {3, 6, 15, 8, 12},
    {4, 12, 7, 15, 9},
    {4, 8, 15, 10, 13},
    {5, 13, 9, 15, 11},
    {5, 10, 15, 6, 14},
    {0, 3, 7, 8, 4, 1},
    {1, 4, 9, 10, 5, 2},
    {0, 2, 5, 11, 6, 3},
    {11, 10, 9, 8, 7, 6}
};

const std::vector<std::vector<unsigned int>>  c30_adjacency = {
    {1, 2, 3, 4, 5},
    {0, 5, 16, 12, 2},
    {0, 1, 12, 13, 3},
    {0, 2, 13, 14, 4},
    {0, 3, 14, 15, 5},
    {0, 4, 15, 16, 1},
    {11, 7, 16, 15, 10},
    {11, 8, 12, 16, 6},
    {11, 9, 13, 12, 7},
    {11, 10, 14, 13, 8},
    {11, 6, 15, 14, 9},
    {10, 9, 8, 7, 6},
    {1, 16, 7, 8, 13, 2},
    {2, 12, 8, 9, 14, 3},
    {3, 13, 9, 10, 15, 4},
    {4, 14, 10, 6, 16, 5},
    {1, 5, 15, 6, 7, 12}
};

inline dual_fullerene create_c20_fullerene() {
    return dual_fullerene(c20_adjacency);
}

inline dual_fullerene create_c28_fullerene() {
    return dual_fullerene(c28_adjacency);
}

inline dual_fullerene create_c30_fullerene() {
    return dual_fullerene(c30_adjacency);
}

#endif //CONSTRUCT_H
