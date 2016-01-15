#include <STL_Example5_Schwarzschild.hpp>

// Simulation params
using integral = std::int32_t;      // Type used to represent index values
using real = float;                 // Type used to represent real values
using complex = std::complex<real>; // Type used to represent complex values
constexpr integral L_max = 7;       // Maximum multipole values for series expansion
constexpr integral S_max = 3;       // Maximum spin values for series expansion

// Convenience type aliases
using namespace Multipole::stl;
using radial_index = integral;
using radial_extent = Radial::Extent<radial_index>;
using spherical_extent = SpinWeightedSpherical::Extent<L_max, S_max, integral>;
using spherical_index = SpinWeightedSpherical::Index<L_max, S_max, integral>;
using gaunt_matrix = SpinWeightedGaunt::Matrix<L_max, S_max, integral, real>;
template <typename T> using coeff_vector = SpinWeightedSpherical::Vector<L_max, S_max, Parity::Even, integral, T>;
template <typename T> using radial_coeff_vector = Radial::Vector<real, integral, coeff_vector<T>>;
using state_vector = PDE::StateVector<radial_coeff_vector<real>, radial_coeff_vector<complex>>;
using solver = PDE::RK4::Solver<radial_coeff_vector<real>, radial_coeff_vector<complex>>;

using time_point = std::chrono::high_resolution_clock::time_point;

int main(int argc, char* argv[])
{
    // Convenience enum for equation numbering
    enum Eq
    {
        K = 0,
        k = 1
    };

    // Declare variables
    std::size_t count;
    real drho, one_per_two;
    radial_extent rho_ext;
    spherical_extent lms_ext;
    radial_coeff_vector<real> N_kalap, K_kalap, kappa;
    radial_coeff_vector<real> K_var;
    radial_coeff_vector<complex> k_var;
    solver rk4;

    // Initialize variables
    count = 1;
    drho = 0.1f;
    rho_ext = radial_extent(0, 1023);
    lms_ext = spherical_extent({ 0, 0, 0 }, { L_max, L_max, S_max });

    rk4.setLHS(state_vector(K_var, k_var));
    rk4.setEquationSystem([&](const state_vector& state)
    {
        return state_vector();
    });

    std::cout << "Done." << std::endl;

    return EXIT_SUCCESS;
}