#include <iostream>
#include <fullerene/construct.h>


int main() {
    auto C20 = create_c20_fullerene();
    std::cout << "Successfully created C20 fullerene" << std::endl;

    auto C28 = create_c28_fullerene();
    std::cout << "Successfully created C28 fullerene" << std::endl;

    auto C30 = create_c30_fullerene();
    std::cout << "Successfully created C30 fullerene" << std::endl;

    auto C20_primal = C20.to_primal();
    std::cout << "Successfully converted C20 dual to primal" << std::endl;

    auto C28_primal = C28.to_primal();
    std::cout << "Successfully converted C28 dual to primal" << std::endl;

    auto C30_primal = C30.to_primal();
    std::cout << "Successfully converted C30 dual to primal" << std::endl;

    return 0;
}
