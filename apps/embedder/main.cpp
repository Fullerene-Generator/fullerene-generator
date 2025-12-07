#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <embeddings/embedder.h>

static graph read_graph_from_stdin() {
    std::size_t n;
    if (!(std::cin >> n)) {
        throw std::runtime_error("Failed to read number of vertices");
    }

    std::array<unsigned, 5> outer{};
    for (int i = 0; i < 5; i++) {
        if (!(std::cin >> outer[i])) {
            throw std::runtime_error("Failed to read outer face nodes");
        }
    }

    std::vector<std::array<unsigned, 3>> adj(n);
    for (std::size_t i = 0; i < n; i++) {
        unsigned a, b, c;
        if (!(std::cin >> a >> b >> c)) {
            throw std::runtime_error("Failed to read adjacency list");
        }
        adj[i] = {a, b, c};
    }

    return graph{adj, outer};
}

static void write_embedding(const std::vector<std::array<double,2>>& coords) {
    for (auto const& p : coords) {
        std::cout << p[0] << " " << p[1] << "\n";
    }
}

int main() {
    try {
        graph g = read_graph_from_stdin();

        auto coords = embedder::compute_tutte(g);

        write_embedding(coords);

    } catch (const std::exception& ex) {
        std::cerr << "Embedder error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
