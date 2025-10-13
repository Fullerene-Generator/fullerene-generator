#ifndef FULLERENE_H
#define FULLERENE_H
#include <memory>
#include <vector>

enum class node_type {
    NODE_5,
    NODE_6,
  };

struct node {
    node_type type;
    std::vector<std::shared_ptr<node>> neighbors;
};

class dual_fullerene {
    std::vector<std::shared_ptr<node>> nodes;

public:
    explicit dual_fullerene(std::vector<std::shared_ptr<node>> nodes);
};

#endif //FULLERENE_H
