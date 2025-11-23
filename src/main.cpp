#include <iostream>
#include <fullerene/construct.h>
#include <fullerene/dual_fullerene.h>
#include <fullerene/dual_fullerene.h>
#include <fullerene/fullerene.h>

#include <expansions/l_expansion.h>



int main() {
    // ----- 1. Create C20 dual -----
    auto C20 = create_c20_fullerene();
    std::cout << "Successfully created C20 fullerene (dual)" << std::endl;

    // ----- 2. Visualize BEFORE any L0 expansion -----
    {
        auto C20_primal = C20.to_primal();
        std::cout << "Converted C20 dual to primal" << std::endl;

        C20_primal.compute_tutte_embedding();
        C20_primal.write_graphviz("C20_before_L0.svg", true);
        std::cout << "Wrote C20_before_L0.svg" << std::endl;
    }

    // ----- 3. Find L0 candidates on the dual -----
    auto candidates = find_L0_candidates(C20);
    std::cout << "Found " << candidates.size() << " L0 candidates" << std::endl;

    if (candidates.empty()) {
        std::cout << "No L0 candidates found on C20, nothing to expand." << std::endl;
        return 0;
    }

    // ----- 4. Take the first candidate and wrap it in an L0_expansion -----
    const auto& cand = candidates[0];
    L0_expansion expansion(C20, cand);

    if (!expansion.validate()) {
        std::cerr << "First L0 candidate failed validate()." << std::endl;
        return 1;
    }

    // Apply expansion directly on the dual C20
    expansion.apply();
    std::cout << "Applied one L0 expansion on C20" << std::endl;

    // ----- 5. Visualize AFTER the L0 expansion -----
    {
        auto C20_after_primal = C20.to_primal();
        std::cout << "Converted expanded dual to primal" << std::endl;

        C20_after_primal.compute_tutte_embedding();
        C20_after_primal.write_graphviz("C20_after_L0.svg", true);
        std::cout << "Wrote C20_after_L0.svg" << std::endl;
    }

    return 0;
}
