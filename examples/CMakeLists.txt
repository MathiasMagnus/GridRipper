###############################################################################
#
# CMakeLists.txt file for the GridRipper example applications
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2016/01/15
#
###############################################################################

# The project name and version, release numbers
project (Examples)
set     (Examples_VERSION_MAJOR 0)
set     (Examples_VERSION_MINOR 1)
set     (Examples_VERSION_PATCH 0)

add_subdirectory (STL-Example1-GauntDistribution)
add_subdirectory (STL-Example2-GauntContraction)
add_subdirectory (STL-Example3-RungeKutta4)
add_subdirectory (STL-Example4-Lattice1D)
add_subdirectory (STL-Example5-Schwarzschild)