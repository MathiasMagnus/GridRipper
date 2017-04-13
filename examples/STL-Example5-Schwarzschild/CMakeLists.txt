###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Example5-Schwarzschild example
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2016/15/01
#
###############################################################################

# The project name and version, release numbers
project (STL-Example5-Schwarzschild)

# Adding source code files according to configuration
set (STL_Example5_Schwarzschild_BUILD)
set (STL_Example5_Schwarzschild_HDRS)
set (STL_Example5_Schwarzschild_SRCS)

list (APPEND STL_Example5_Schwarzschild_HDRS ${PROJECT_SOURCE_DIR}/inc/STL_Example5_Schwarzschild.hpp)

list (APPEND STL_Example5_Schwarzschild_SRCS ${PROJECT_SOURCE_DIR}/src/STL_Example5_Schwarzschild.cpp)

set (STL_Example5_Schwarzschild_BUILD ${STL_Example5_Schwarzschild_HDRS}
                                      ${STL_Example5_Schwarzschild_SRCS})

# Adding library target for build
add_executable (${PROJECT_NAME} ${STL_Example5_Schwarzschild_BUILD})

# Adding extra library dependencies
target_link_libraries (${PROJECT_NAME} ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc/
                                                   ${Gripper_ADDITIONAL_INCLUDE_DIRS})

# Create project groups for IDEs
set_target_properties (${PROJECT_NAME} PROPERTIES FOLDER "Examples//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Example5_Schwarzschild_HDRS})
source_group ("Sources" FILES ${STL_Example5_Schwarzschild_SRCS})