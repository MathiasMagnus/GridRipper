#include <STL-Test3-RK4.hpp>


int main(int argc, char* argv[])
{
    // Initialize Gripper runtime
    Gripper::stl::initialize(argc, argv);

    // Start time
    auto start = std::chrono::high_resolution_clock::now();

    // Utility typedef to save some typing
    typedef PDE::StateVector<double, float> States;

    // Create a solver for a given set of variables
    PDE::RK4::Solver<double, float> solver;

    // Set the equation system to integrate
    //
    // The equation system must be a callable that has an input of inmutable set of states and returns with a new set of states
    solver.setEquationSystem([](const States& states) -> States
    {
        return States(std::exp(states.get<0>()), std::exp(states.get<1>()));
    });

    // Set the initial values for the solver
    //
    // After providing the solver with the initial conditions, it will store it internally, hence setLHS requires RValue reference
    solver.setLHS(States(static_cast<double>(M_E), static_cast<float>(M_E)));

    // Do the actual time integration
    float dt = 0.01f;
    for(float t = 0.0f; t < 1.0f; t += dt) solver.iterate(dt);

    // Check if results are valid
    int result =
        (std::abs(solver.getLHS().get<0>() - static_cast<double>(M_E)) < 1e-8) &&
        (std::abs(solver.getLHS().get<1>() - static_cast<float>(M_E)) < 1e-6)
        ? 0 : -1;
        
    // Stop time
    auto end = std::chrono::high_resolution_clock::now();

    Gripper::stl::clog << Gripper::TIMING << "Tests took " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds.";

    if(result)
        Gripper::stl::clog << Gripper::INFO << "Test succeeded.";
    else
        Gripper::stl::clog << Gripper::INFO << "Test failed.";
	
	return result;
}