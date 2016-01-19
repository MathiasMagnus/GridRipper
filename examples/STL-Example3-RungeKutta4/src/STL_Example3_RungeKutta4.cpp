#include <STL_Example3_RungeKutta4.hpp>

#include <complex>

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
    using state_vector = PDE::StateVector<real, coeff_vector, coeff_vector>;
    using solver = PDE::RK4::Solver<real, state_vector>;

    // Declare variables
    time_point start, end;
    extent ext(index{ 0, 0, 0 }, index{ L_max, L_max, S_max });
    index ai{ 3, -3, 1 }, bi{ 3, 3, -1 }, ci{ 0, 0, 0 };
    coeff_vector c(ext);
    gaunt_matrix gaunt;
    real a, b, dt;
    solver rk4;
    std::size_t count;

    // Initialize Gaunt-matrix
    std::cout << "Calculating gaunt spin-weighted Gaunt-matrix of L_max = " << static_cast<int>(L_max) << " and S_max = " << static_cast<int>(S_max) << std::endl;
    start = std::chrono::high_resolution_clock::now();

    gaunt.calculate();

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Count of non-zero elements = " << gaunt.size() << std::endl;
    std::cout << "Calculation took " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " sec\n" << std::endl;

    // Initialize state vector
    for (auto& coeff : c) coeff = 1;
    a = 1;
    b = 0.01f;
    dt = 0.1f;
    count = 100;

    rk4.setLHS(state_vector(a, c, c));
    rk4.setEquationSystem([&](const state_vector& state)
    {
        return state_vector(
            std::exp(-b * state.get<0>()),
            std::exp(-b * SpinWeightedGaunt::contract(gaunt, state.get<1>(), state.get<1>())),
            eth_bar(eth(state.get<2>())));
    });

    // Integrate given equation
    std::cout << "Performing " << count << " RK4 steps with dt = " << dt << std::endl;
    start = std::chrono::high_resolution_clock::now();

    for (std::size_t n = 0; n < count; ++n)
        rk4.iterate(dt);

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Integration took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n" << std::endl;

    return EXIT_SUCCESS;
}
