#include <STL_Example4_Lattice1D.hpp>

int main(int argc, char* argv[])
{
    // Simulation params
    using integral = std::int32_t;  // Type used to represent index values
    using real = float;             // Type used to represent real values

    // Convenience type aliases
    using namespace Multipole::stl;
    using index = integral;
    using extent = Radial::Extent<integral>;
    using vector = Radial::Vector<real, index, real>;

    // Declare variables
    real dr;
    extent ext;
    vector vec;

    // Initialize variables
    dr = 0.1f;
    ext = extent(0, 2047);
    vec = vector(dr, ext);
    auto r_sq = Radial::func(dr, ext, [](const index& i) { return static_cast<real>(i*i); });

    // Calculate a simple equation
    vec = derivate(-r_sq);

    std::cout << "Done." << std::endl;

    return EXIT_SUCCESS;
}
