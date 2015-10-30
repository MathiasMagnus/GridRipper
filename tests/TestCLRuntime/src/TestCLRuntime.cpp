#include <TestSTLRuntime.hpp>

int main(int argc, char* argv[])
{
    Gripper::stl::initialize(argc, argv);

    if (Gripper::stl::gaunt.size() == 0) return EXIT_FAILURE;
	
	return EXIT_SUCCESS;
}