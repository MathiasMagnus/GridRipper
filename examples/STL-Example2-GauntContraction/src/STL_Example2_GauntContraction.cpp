#include <STL_Example2_GauntContraction.hpp>

int main(int argc, char* argv[])
{
    // Simulation params
    using integral = std::int32_t;  // Type used to represent index values
    using real = float;             // Type used to represent real values
    constexpr integral L_max = 7;   // Maximum multipole values for series expansion
    constexpr integral S_max = 3;   // Maximum spin values for series expansion

    // Convenience type aliases
    using namespace Multipole::stl;
    using extent = SpinWeightedSpherical::Extent<L_max, S_max, integral>;
    using index = SpinWeightedSpherical::Index<L_max, S_max, integral>;
    using gaunt_matrix = SpinWeightedGaunt::Matrix<L_max, S_max, integral, real>;
    using coeff_vector = SpinWeightedSpherical::Vector<L_max, S_max, Parity::Even, integral, real>;
    using time_point = std::chrono::high_resolution_clock::time_point;

    // Declare variables
    time_point start, end;
    extent ext(index{0, 0, 0}, index{L_max, L_max, S_max});
    index ai{ 3, -3, 1 }, bi{ 3, 3, -1 }, ci{ 0, 0, 0 };
    coeff_vector a(ext), b(ext), c;
    gaunt_matrix gaunt;

    // Initialize Gaunt-matrix
    std::cout << "Calculating gaunt spin-weighted Gaunt-matrix of L_max = " << static_cast<int>(L_max) << " and S_max = " << static_cast<int>(S_max) << std::endl;
    start = std::chrono::high_resolution_clock::now();

    gaunt.calculate();

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Count of non-zero elements = " << gaunt.size() << std::endl;
    std::cout << "Calculation took " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " sec\n" << std::endl;

    // Initialize coefficient vectors
    a.at(ai) = 1;
    b.at(bi) = 1;
        
    // Contract Gaunt-matrix with an expansion of a field
    std::cout << "Contracting gaunt matrix with a" << ai << " = " << 1 << " and b" << bi << " = " << 1 << std::endl;
    start = std::chrono::high_resolution_clock::now();

    c = SpinWeightedGaunt::contract(gaunt, a, b);

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Result at c" << ci << " = " << c.at(ci) << std::endl;
    std::cout << "Contraction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n" << std::endl;
    
	return EXIT_SUCCESS;
}
