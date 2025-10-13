#include "fullerene.h"

#include <utility>

dual_fullerene::dual_fullerene(std::vector<std::shared_ptr<node>> nodes) {
    this->nodes = std::move(nodes);
}
