###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test3-RK4 test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 23/07/2014
#
###############################################################################

# The project name and version, release numbers
project (STL-Test3-RK4)


# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
include_directories (${PROJECT_SOURCE_DIR}/inc/)
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (STL_Test3_RK4_BUILD)
set (STL_Test3_RK4_HDRS)
set (STL_Test3_RK4_SRCS)

set (STL_Test3_RK4_HDRS ${STL_Test3_RK4_HDRS} ${PROJECT_SOURCE_DIR}/inc/STL-Test3-RK4.hpp)

set (STL_Test3_RK4_SRCS ${STL_Test3_RK4_SRCS} ${PROJECT_SOURCE_DIR}/src/STL-Test3-RK4.cpp)

set (STL_Test3_RK4_BUILD ${STL_Test3_RK4_BUILD} ${STL_Test3_RK4_HDRS} ${STL_Test3_RK4_SRCS})


# Build and link
add_executable( STL-Test3-RK4 ${STL_Test3_RK4_BUILD})
target_link_libraries( STL-Test3-RK4 Gripper ${Gripper_ADDITIONAL_LIBS})


# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test3_RK4_HDRS})
source_group ("Sources" FILES ${STL_Test3_RK4_SRCS})

# Create project groups for IDEs
set_target_properties (STL-Test3-RK4 PROPERTIES FOLDER "Tests//STL")


# Add CTest entry
add_test( STL-Test3-RK4 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/STL-Test3-RK4 )


# Platform-specific compiler options
if(MSVC_VERSION)
	set_target_properties (STL-Test3-RK4 PROPERTIES COMPILE_FLAGS ${MSVC_COMPILER_FLAGS})
	set_target_properties (STL-Test3-RK4 PROPERTIES LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)