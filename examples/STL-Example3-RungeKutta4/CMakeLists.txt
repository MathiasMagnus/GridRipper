###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Example3-RungeKutta4 example
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2016/05/01
#
###############################################################################

# The project name and version, release numbers
project (STL-Example3-RungeKutta4)


# Adding source code files according to configuration
set (STL_Example3_RungeKutta4_BUILD)
set (STL_Example3_RungeKutta4_HDRS)
set (STL_Example3_RungeKutta4_SRCS)

list (APPEND STL_Example3_RungeKutta4_HDRS ${PROJECT_SOURCE_DIR}/inc/STL_Example3_RungeKutta4.hpp)

list (APPEND STL_Example3_RungeKutta4_SRCS ${PROJECT_SOURCE_DIR}/src/STL_Example3_RungeKutta4.cpp)

set (STL_Example3_RungeKutta4_BUILD ${STL_Example3_RungeKutta4_HDRS}
                                          ${STL_Example3_RungeKutta4_SRCS})


# Adding library target for build
add_executable (STL-Example3-RungeKutta4 ${STL_Example3_RungeKutta4_BUILD})

# Adding extra library dependencies
target_link_libraries (STL-Example3-RungeKutta4 ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc/
                                                   ${Gripper_ADDITIONAL_INCLUDE_DIRS})

# Create project groups for IDEs
set_target_properties (STL-Example3-RungeKutta4 PROPERTIES FOLDER "Examples//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Example3_RungeKutta4_HDRS})
source_group ("Sources" FILES ${STL_Example3_RungeKutta4_SRCS})