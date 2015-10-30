#include <STL-Test4-Plotter.hpp>


int main(int argc, char* argv[])
{
    // Using directives
    using namespace Multipole::stl;

    // Initialize Gripper runtime
    Gripper::stl::initialize(argc, argv);

    // Start time
    auto start = std::chrono::high_resolution_clock::now();

    // Allocate memory
    Radial::Vector<float> a(Radial::Extent(0, 1024));
    Radial::Vector<float> b(Radial::Extent(0, 1024));

    // Initial data
    for (auto r = a.extent().initial(); a.extent().contains(r); ++r)
    {
        a.at(r) = 1.0f * std::exp(-std::pow(static_cast<float>((r - (a.extent().final() - a.extent().initial())/2).r), 2) / (2 * std::pow(64.f, 2)));
        b.at(r) = 0.5f * std::exp(-std::pow(static_cast<float>((r - (b.extent().final() - b.extent().initial())/2).r), 2) / (2 * std::pow(128.f, 2)));
    }

    Gripper::stl::GnuPlotter plotter("plot");

    plotter << a << b;

    plotter.close();
        
    // Stop time
    auto end = std::chrono::high_resolution_clock::now();

    Gripper::stl::clog << Gripper::TIMING << "Tests took " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds.";
	
	return 0;
}