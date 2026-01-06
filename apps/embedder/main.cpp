#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <embeddings/embedder.h>

#include "fullerene/construct.h"

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

static void write_embedding_2d(const std::vector<std::array<double,2>>& coords) {
    for (auto const& p : coords) {
        std::cout << p[0] << " " << p[1] << "\n";
    }
}

static void write_embedding_3d(const std::vector<std::array<double,3>>& coords) {
    for (auto const& p : coords) {
        std::cout << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "Usage: " << argv[0] << " <2|3>\n";
            return 1;
        }

        int mode = std::stoi(argv[1]);
        if (mode != 2 && mode != 3) {
            std::cerr << "Mode must be 2 or 3\n";
            return 1;
        }

        graph g = read_graph_from_stdin();

        if (mode == 2) {
            const auto coords = embedder::compute_tutte(g);
            write_embedding_2d(coords);
        } else {
            const auto coords = embedder::compute_3d_force_embedding(g);
            write_embedding_3d(coords);
        }
    } catch (const std::exception& ex) {
        std::cerr << "Embedder error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
