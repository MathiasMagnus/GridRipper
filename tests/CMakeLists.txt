###############################################################################
#
# CMakeLists.txt file for the GridRipper automated tests
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 23/17/2014
#
###############################################################################

# The project name and version, release numbers
project (Tests)
set     (Tests_VERSION_MAJOR 0)
set     (Tests_VERSION_MINOR 5)
set     (Tests_VERSION_PATCH 0)

# Process tests
add_subdirectory (Base-Test1-Runtime)

if (USE_STL)
	add_subdirectory (STL-Test1-Runtime)
	add_subdirectory (STL-Test2-Arithmetics)
	add_subdirectory (STL-Test3-RK4)
	add_subdirectory (STL-Test4-Plotter)
endif (USE_STL)

if (USE_OPENCL)
	add_subdirectory (TestCLRuntime)
	add_subdirectory (TestCLArithmetics)
endif (USE_OPENCL)