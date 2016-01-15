###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Example2-GauntContraction example
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2016/05/01
#
###############################################################################

# The project name and version, release numbers
project (STL-Example2-GauntContraction)


# Adding source code files according to configuration
set (STL_Example2_GauntContraction_BUILD)
set (STL_Example2_GauntContraction_HDRS)
set (STL_Example2_GauntContraction_SRCS)

list (APPEND STL_Example2_GauntContraction_HDRS ${PROJECT_SOURCE_DIR}/inc/STL_Example2_GauntContraction.hpp)

list (APPEND STL_Example2_GauntContraction_SRCS ${PROJECT_SOURCE_DIR}/src/STL_Example2_GauntContraction.cpp)

set (STL_Example2_GauntContraction_BUILD ${STL_Example2_GauntContraction_HDRS}
                                          ${STL_Example2_GauntContraction_SRCS})


# Adding library target for build
add_executable (STL-Example2-GauntContraction ${STL_Example2_GauntContraction_BUILD})

# Adding extra library dependencies
target_link_libraries (STL-Example2-GauntContraction ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc/
                                                   ${Gripper_ADDITIONAL_INCLUDE_DIRS})

# Create project groups for IDEs
set_target_properties (STL-Example2-GauntContraction PROPERTIES FOLDER "Examples//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Example2_GauntContraction_HDRS})
source_group ("Sources" FILES ${STL_Example2_GauntContraction_SRCS})