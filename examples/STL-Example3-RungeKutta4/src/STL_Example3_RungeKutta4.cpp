#include <STL_Example3_RungeKutta4.hpp>

enum Population : std::size_t
{
    X = 0,
    Y = 1
};

int main(int argc, char* argv[])
{
    // Simulation params
    using integral = std::int32_t;  // Type used to represent index values
    using real = float;             // Type used to represent real values

    // Convenience type aliases
    using time_point = std::chrono::high_resolution_clock::time_point;
    using state_vector = PDE::StateVector<real, real>;
    using solver = PDE::RK4::Solver<real, state_vector>;

    // Declare variables
    time_point start, end;
    real a, b, c, d, x, y, dt;
    solver rk4;
    std::size_t count;

    // Initialize simulation params
    x = 7;          // Rabbit population
    y = 5;          // Fox population

    a = 0.35f;      // Rabbit multiplicative factor
    b = 0.15f;      // De-rabitization factor
    c = 0.05f;      // Fox growth factor
    d = 0.15f;      // Fox aging

    dt = 1.0f;      // Time step
    count = 500;    // Number of steps to take

    rk4.lhs() = state_vector(x, y);
    rk4.equation() = [&](state_vector& result, const state_vector& state)
    {
        result = PDE::make_equation(a * state.get<X>()                  - b * state.get<X>() * state.get<Y>(),
                                    c * state.get<X>() * state.get<Y>() - d * state.get<Y>());
    };

    // Integrate given equation
    std::cout << "Performing " << count << " RK4 steps with dt = " << dt << std::endl;
    start = std::chrono::high_resolution_clock::now();

    std::ostringstream buffer;

    for (std::size_t n = 0; n < count; ++n)
    {
        rk4.iterate(dt);
        buffer << n * dt << " " << rk4.lhs().get<X>() << " " << rk4.lhs().get<Y>() << std::endl;
    }

    std::ofstream output("lotka_volterra.dat");

    output << buffer.str();

    output.close();

    end = std::chrono::high_resolution_clock::now();
    std::cout << "Integration took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n" << std::endl;

    return EXIT_SUCCESS;
}
