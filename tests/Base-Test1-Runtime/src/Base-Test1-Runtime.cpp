#include <Base-Test1-Runtime.hpp>

int main(int argc, char* argv[])
{
    Gripper::initialize(argc, argv);

    if (Gripper::gaunt.size() == 0) return EXIT_FAILURE;
	
	return EXIT_SUCCESS;
}