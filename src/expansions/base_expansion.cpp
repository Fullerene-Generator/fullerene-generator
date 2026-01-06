#include <expansions/base_expansion.h>
#include <unordered_set>
#include <vector>

bool patch_nodes_unique(const std::vector<int>& path, const std::vector<int>& parallel_path) {
    std::unordered_set<int> seen;
    seen.reserve(path.size() + parallel_path.size());

    for (int v : path) {
        if (!seen.insert(v).second) {
            return false;
        }
    }
    for (int v : parallel_path) {
        if (!seen.insert(v).second) {
            return false;
        }
    }
    return true;
}
