#include <STL_Example1_GauntDistribution.hpp>

#include <cstdint>
#include <iostream>

int main(int argc, char* argv[])
{
	Gripper::stl::initialize(argc, argv);

    // Convenience type aliases
    using namespace Multipole::stl;
    using index_type = std::int32_t;
    using floating_type = float;
    using extent_type = SpinWeightedSpherical::Extent<index_type>;
    using gaunt_matrix = SpinWeightedGaunt::Matrix<index_type, floating_type>;
    using spherical_vector = SpinWeightedSpherical::Vector<index_type, floating_type>;

    // Simulation params
    constexpr index_type max_L = 7;
    constexpr index_type max_S = 3;

    // Initialize Gaunt-matrix
    gaunt_matrix gaunt(max_L, max_S);

    std::cout << std::endl;
    std::cout << "L_max = " << max_L << "\nS_max = " << max_S << "\n";
    std::cout << "Size of spin-weighted gaunt matrix = " << gaunt.size() << std::endl;

    spherical_vector a(max_L, max_S);
    spherical_vector b(max_L, max_S);

    a.at(spherical_vector::index_type{ 3, -3, 1 }) = 1;
    b.at(spherical_vector::index_type{ 3, 3, -1 }) = 1;

    std::cout << a.at(spherical_vector::index_type{ 3, -3, 1 }) << std::endl;
    std::cout << b.at(spherical_vector::index_type{ 3, 3, -1 }) << std::endl;
        
    // Contract Gaunt-matrix with an expansion of a field
    auto start = std::chrono::high_resolution_clock::now();
    spherical_vector c = SpinWeightedGaunt::contract(gaunt, a, b);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << c.at(spherical_vector::index_type{ 0, 0, 0 }) << std::endl;
    std::cout << "Contraction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    
	return EXIT_SUCCESS;
}