###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test2-Arithmetics test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 23/07/2014
#
###############################################################################

# The project name and version, release numbers
project (STL-Test2-Arithmetics)


# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
include_directories (${PROJECT_SOURCE_DIR}/inc/)
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (STL_Test2_Arithmetics_BUILD)
set (STL_Test2_Arithmetics_HDRS)
set (STL_Test2_Arithmetics_SRCS)

set (STL_Test2_Arithmetics_HDRS ${STL_Test2_Arithmetics_HDRS} ${PROJECT_SOURCE_DIR}/inc/STL-Test2-Arithmetics.hpp)

set (STL_Test2_Arithmetics_SRCS ${STL_Test2_Arithmetics_SRCS} ${PROJECT_SOURCE_DIR}/src/STL-Test2-Arithmetics.cpp)

set (STL_Test2_Arithmetics_BUILD ${STL_Test2_Arithmetics_BUILD} ${STL_Test2_Arithmetics_HDRS} ${STL_Test2_Arithmetics_SRCS})


# Build and link
add_executable( STL-Test2-Arithmetics ${STL_Test2_Arithmetics_BUILD})
target_link_libraries( STL-Test2-Arithmetics Gripper ${Gripper_ADDITIONAL_LIBS})


# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test2_Arithmetics_HDRS})
source_group ("Sources" FILES ${STL_Test2_Arithmetics_SRCS})

# Create project groups for IDEs
set_target_properties (STL-Test2-Arithmetics PROPERTIES FOLDER "Tests//STL")


# Add CTest entry
add_test( STL-Test2-Arithmetics ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/STL-Test2-Arithmetics )


# Platform-specific compiler options
if(MSVC_VERSION)
	set_target_properties (STL-Test2-Arithmetics PROPERTIES COMPILE_FLAGS ${MSVC_COMPILER_FLAGS})
	set_target_properties (STL-Test2-Arithmetics PROPERTIES LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)