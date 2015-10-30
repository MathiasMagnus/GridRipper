###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test4-Plotter test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 23/07/2014
#
###############################################################################

# The project name and version, release numbers
project (STL-Test4-Plotter)


# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
include_directories (${PROJECT_SOURCE_DIR}/inc/)
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (STL_Test4_Plotter_BUILD)
set (STL_Test4_Plotter_HDRS)
set (STL_Test4_Plotter_SRCS)

set (STL_Test4_Plotter_HDRS ${STL_Test4_Plotter_HDRS} ${PROJECT_SOURCE_DIR}/inc/STL-Test4-Plotter.hpp)

set (STL_Test4_Plotter_SRCS ${STL_Test4_Plotter_SRCS} ${PROJECT_SOURCE_DIR}/src/STL-Test4-Plotter.cpp)

set (STL_Test4_Plotter_BUILD ${STL_Test4_Plotter_BUILD} ${STL_Test4_Plotter_HDRS} ${STL_Test4_Plotter_SRCS})


# Build and link
add_executable( STL-Test4-Plotter ${STL_Test4_Plotter_BUILD})
target_link_libraries( STL-Test4-Plotter Gripper ${Gripper_ADDITIONAL_LIBS})


# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test4_Plotter_HDRS})
source_group ("Sources" FILES ${STL_Test4_Plotter_SRCS})

# Create project groups for IDEs
set_target_properties (STL-Test4-Plotter PROPERTIES FOLDER "Tests//STL")


# Add CTest entry
add_test( STL-Test4-Plotter ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/STL-Test4-Plotter )


# Platform-specific compiler options
if(MSVC_VERSION)
	set_target_properties (STL-Test4-Plotter PROPERTIES COMPILE_FLAGS ${MSVC_COMPILER_FLAGS})
	set_target_properties (STL-Test4-Plotter PROPERTIES LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)