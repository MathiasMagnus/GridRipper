###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test5-Arithmetics-SW-Compound test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 04/02/2016
#
###############################################################################

# The project name and version, release numbers
project (STL-Test5-Arithmetics-SW-Compound)


# Adding source code files according to configuration
set (STL_Test5_Arithmetics_SW_Compound_BUILD)
set (STL_Test5_Arithmetics_SW_Compound_HDRS)
set (STL_Test5_Arithmetics_SW_Compound_SRCS)

list (APPEND STL_Test5_Arithmetics_SW_Compound_HDRS ${PROJECT_SOURCE_DIR}/inc/STL_Test5_Arithmetics_SW_Compound.hpp)

list (APPEND STL_Test5_Arithmetics_SW_Compound_SRCS ${PROJECT_SOURCE_DIR}/src/STL_Test5_Arithmetics_SW_Compound.cpp)

set (STL_Test5_Arithmetics_SW_Compound_BUILD ${STL_Test5_Arithmetics_SW_Compound_HDRS}
                                          ${STL_Test5_Arithmetics_SW_Compound_SRCS})


# Adding library target for build
add_executable (STL-Test5-Arithmetics-SW-Compound ${STL_Test5_Arithmetics_SW_Compound_BUILD})

# Adding extra library dependencies
target_link_libraries (${PROJECT_NAME} ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc/
                                                   ${Gripper_ADDITIONAL_INCLUDE_DIRS})

# Create project groups for IDEs
set_target_properties (STL-Test5-Arithmetics-SW-Compound PROPERTIES FOLDER "Tests//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test5_Arithmetics-SW-Compound_HDRS})
source_group ("Sources" FILES ${STL_Test5_Arithmetics-SW-Compound_SRCS})


# Add CTest entry
add_test( ${PROJECT_NAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/STL-Test5-Arithmetics-SW-Compound )