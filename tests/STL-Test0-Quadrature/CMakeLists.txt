###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test1-RK4 test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 04/02/2016
#
###############################################################################

# The project name and version, release numbers
project (STL-Test0-Quadrature)

# Adding source code files according to configuration
set (STL_Test0_Quadrature_BUILD)
set (STL_Test0_Quadrature_HDRS)
set (STL_Test0_Quadrature_SRCS)

list (APPEND STL_Test0_Quadrature_HDRS ${PROJECT_SOURCE_DIR}/inc/STL_Test0_Quadrature.hpp)

list (APPEND STL_Test0_Quadrature_SRCS ${PROJECT_SOURCE_DIR}/src/STL_Test0_Quadrature.cpp)

set (STL_Test0_Quadrature_BUILD ${STL_Test0_Quadrature_HDRS}
                                ${STL_Test0_Quadrature_SRCS})

# Adding library target for build
add_executable (${PROJECT_NAME} ${STL_Test0_Quadrature_BUILD})

# Adding extra library dependencies
target_link_libraries (${PROJECT_NAME} ${Gripper_ADDITIONAL_LIBS})

# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
target_include_directories (${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/inc/
                                                   ${Gripper_ADDITIONAL_INCLUDE_DIRS})
# Enable C++11 language support
set_property (TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)
set_property (TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

# Create project groups for IDEs
set_target_properties (${PROJECT_NAME} PROPERTIES FOLDER "Tests//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test0_Quadrature_HDRS})
source_group ("Sources" FILES ${STL_Test0_Quadrature_SRCS})

# Add CTest entry
add_test (NAME ${PROJECT_NAME}
          COMMAND ${PROJECT_NAME})