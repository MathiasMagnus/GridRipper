###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Example4-Lattice1D example
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2016/05/01
#
###############################################################################

# The project name and version, release numbers
project (STL-Example4-Lattice1D)


# Adding source code files according to configuration
set (STL_Example4_Lattice1D_BUILD)
set (STL_Example4_Lattice1D_HDRS)
set (STL_Example4_Lattice1D_SRCS)

list (APPEND STL_Example4_Lattice1D_HDRS ${PROJECT_SOURCE_DIR}/inc/STL_Example4_Lattice1D.hpp)

list (APPEND STL_Example4_Lattice1D_SRCS ${PROJECT_SOURCE_DIR}/src/STL_Example4_Lattice1D.cpp)

set (STL_Example4_Lattice1D_BUILD ${STL_Example4_Lattice1D_HDRS}
                                          ${STL_Example4_Lattice1D_SRCS})


# Adding library target for build
add_executable (STL-Example4-Lattice1D ${STL_Example4_Lattice1D_BUILD})

# Adding extra library dependencies
target_link_libraries (STL-Example4-Lattice1D ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc/
                                                   ${Gripper_ADDITIONAL_INCLUDE_DIRS})

# Create project groups for IDEs
set_target_properties (STL-Example4-Lattice1D PROPERTIES FOLDER "Examples//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Example4_Lattice1D_HDRS})
source_group ("Sources" FILES ${STL_Example4_Lattice1D_SRCS})