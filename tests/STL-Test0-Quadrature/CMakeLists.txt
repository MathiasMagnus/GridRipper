###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test1-RK4 test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 04/02/2016
#
###############################################################################

# The name of the target
set (Target_Name STL-Test0-Quadrature)

# Adding source code files according to configuration

set (STL_Test0_Quadrature_HDRS inc/STL_Test0_Quadrature.hpp)

set (STL_Test0_Quadrature_SRCS src/STL_Test0_Quadrature.cpp)

set (STL_Test0_Quadrature_BUILD ${STL_Test0_Quadrature_HDRS}
                                ${STL_Test0_Quadrature_SRCS})

# Adding library target for build
add_executable (${Target_Name} ${STL_Test0_Quadrature_BUILD})

# Adding extra library dependencies
target_link_libraries (${Target_Name} ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${Target_Name} PUBLIC inc/
                                                  ${Gripper_ADDITIONAL_INCLUDE_DIRS})
# Enable C++14 language support
#        project groups for IDEs
set_target_properties (${Target_Name} PROPERTIES CXX_STANDARD 14
                                                 CXX_STANDARD_REQUIRED ON
											     FOLDER "Tests//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test0_Quadrature_HDRS})
source_group ("Sources" FILES ${STL_Test0_Quadrature_SRCS})

# Add CTest entry
add_test (NAME ${Target_Name}
          COMMAND ${Target_Name})