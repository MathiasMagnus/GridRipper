###############################################################################
#
# CMakeLists.txt file for the GridRipper automated tests
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 23/17/2014
#
###############################################################################

# The project name and version, release numbers
#project (Tests)
#set (Tests_VERSION_MAJOR 0)
#set (Tests_VERSION_MINOR 5)
#set (Tests_VERSION_PATCH 0)

# Process tests

if (GR_USE_STL)
  add_subdirectory (STL-Test0-Quadrature)
  #add_subdirectory (STL-Test1-RK4)
  #add_subdirectory (STL-Test2-Arithmetics-Simple)
  #add_subdirectory (STL-Test3-Arithmetics-Compound)
  #add_subdirectory (STL-Test4-Arithmetics-SW-Simple)
  #add_subdirectory (STL-Test5-Arithmetics-SW-Compound)
endif (GR_USE_STL)