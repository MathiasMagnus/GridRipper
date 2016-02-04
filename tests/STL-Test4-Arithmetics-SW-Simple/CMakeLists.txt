###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test2-Arithmetics-SW-Simple test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 04/02/2016
#
###############################################################################

# The project name and version, release numbers
project (STL-Test2-Arithmetics-SW-Simple)


# Adding source code files according to configuration
set (STL_Test2_Arithmetics_SW_Simple_BUILD)
set (STL_Test2_Arithmetics_SW_Simple_HDRS)
set (STL_Test2_Arithmetics_SW_Simple_SRCS)

list (APPEND STL_Test2_Arithmetics_SW_Simple_HDRS ${PROJECT_SOURCE_DIR}/inc/STL_Test2_Arithmetics_SW_Simple.hpp)

list (APPEND STL_Test2_Arithmetics_SW_Simple_SRCS ${PROJECT_SOURCE_DIR}/src/STL_Test2_Arithmetics_SW_Simple.cpp)

set (STL_Test2_Arithmetics_SW_Simple_BUILD ${STL_Test2_Arithmetics_SW_Simple_HDRS}
                                           ${STL_Test2_Arithmetics_SW_Simple_SRCS})


# Adding library target for build
add_executable (STL-Test2-Arithmetics-SW-Simple ${STL_Test2_Arithmetics_SW_Simple_BUILD})

# Adding extra library dependencies
target_link_libraries (${PROJECT_NAME} ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc/
                                                   ${Gripper_ADDITIONAL_INCLUDE_DIRS})

# Create project groups for IDEs
set_target_properties (STL-Test2-Arithmetics-SW-Simple PROPERTIES FOLDER "Tests//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test2_Arithmetics-SW-Simple_HDRS})
source_group ("Sources" FILES ${STL_Test2_Arithmetics-SW-Simple_SRCS})


# Add CTest entry
add_test( ${PROJECT_NAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/STL-Test2-Arithmetics-SW-Simple )