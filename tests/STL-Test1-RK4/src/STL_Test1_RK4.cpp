#include <STL_Test1_RK4.hpp>

int main()
{
    using integral = std::int32_t;          // Type used to represent integral values
    using real = double;                    // Type used to represent real values
    using complex = std::complex<real>;     // Type used to represent complex values
    using vector = std::valarray<real>;     // Type used to represent vector values

    using state_vector = PDE::StateVector<real, complex, vector>;       // Type used as a state vector
    using solver = PDE::RK4::Solver<real, state_vector>;                // Type used to integrate a PDE over state vectors

    using time_point = std::chrono::high_resolution_clock::time_point;  // Time point as returned by times

    // Declare variables    
    real a, x, phi, dx, tolerance;
    complex b;
    vector c;

    // Initialize variables
    a = 100.f;
    b = { 1.f, 0.f };
    c = { 10.f, 50.f, 100.f, 500.f, 1000.f };
    phi = 3.14f / 180;

    x = 4.f;
    dx = 0.1f;
    tolerance = 1e-4f;

    // Initialize solver
    solver rk4;
    rk4.setLHS(state_vector(a, b, c));
    rk4.setEquationSystem([=](const state_vector& state)
    {
        using namespace std::literals;

        return state_vector(
            -state.get<0>(),
            std::exp(complex(0, 1) * phi),
            -state.get<2>()
            );
    });

    for (real X = 0; X < x; X += dx)
    {
        rk4.iterate(dx);
    }

    if (std::abs(rk4.getLHS().get<0>() - a * std::exp(-x)) > tolerance) return EXIT_FAILURE;
    if (false /*TODO: Come up derivate of complex rotation*/) return EXIT_FAILURE;
    if (std::abs(rk4.getLHS().get<2>() - c * std::exp(-x)).max() > tolerance) return EXIT_FAILURE;

	return EXIT_SUCCESS;
}