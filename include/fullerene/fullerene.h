#ifndef FULLERENE_H
#define FULLERENE_H
#include <array>
#include <vector>

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class fullerene {
    std::vector<std::array<unsigned int, 3>> adjacency_;
    std::array<unsigned int, 5> outer_face_nodes_;
    std::vector<std::array<double, 2>> embedding_2d_;

public:
    explicit fullerene(const std::vector<std::array<unsigned int, 3>>& adjacency,
                        const std::array<unsigned int, 5> &outer_face):
                        adjacency_(adjacency),
                        outer_face_nodes_(outer_face) {};

    void compute_tutte_embedding();

    [[nodiscard]] bool has_2d_embedding() const noexcept { return !embedding_2d_.empty(); }
    [[nodiscard]] std::string write_data() const noexcept;
    friend std::ostream &operator<<(std::ostream &os, const fullerene &f);
};

#endif //FULLERENE_H
