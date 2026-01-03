#include <expansions/base_expansion.h>
#include <unordered_set>
#include <vector>

bool patch_nodes_unique(const std::vector<int>& path, const std::vector<int>& para) {
    std::unordered_set<int> seen;
    seen.reserve(path.size() + para.size());

    for (int v : path) {
        if (!seen.insert(v).second) {
            return false;
        }
    }
    for (int v : para) {
        if (!seen.insert(v).second) {
            return false;
        }
    }
    return true;
}
