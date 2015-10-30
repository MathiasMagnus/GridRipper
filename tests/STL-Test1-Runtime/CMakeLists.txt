###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Test1-Runtime test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 23/07/2014
#
###############################################################################

# The project name and version, release numbers
project (STL-Test1-Runtime)


# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
include_directories (${PROJECT_SOURCE_DIR}/inc/)
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (STL_Test1_Runtime_BUILD)
set (STL_Test1_Runtime_HDRS)
set (STL_Test1_Runtime_SRCS)

set (STL_Test1_Runtime_HDRS ${STL_Test1_Runtime_HDRS} ${PROJECT_SOURCE_DIR}/inc/STL-Test1-Runtime.hpp)

set (STL_Test1_Runtime_SRCS ${STL_Test1_Runtime_SRCS} ${PROJECT_SOURCE_DIR}/src/STL-Test1-Runtime.cpp)

set (STL_Test1_Runtime_BUILD ${STL_Test1_Runtime_BUILD} ${STL_Test1_Runtime_HDRS} ${STL_Test1_Runtime_SRCS})


# Build and link
add_executable( STL-Test1-Runtime ${STL_Test1_Runtime_BUILD})
target_link_libraries( STL-Test1-Runtime Gripper ${Gripper_ADDITIONAL_LIBS})


# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Test1_Runtime_HDRS})
source_group ("Sources" FILES ${STL_Test1_Runtime_SRCS})

# Create project groups for IDEs
set_target_properties (STL-Test1-Runtime PROPERTIES FOLDER "Tests//STL")


# Add CTest entry
add_test( STL-Test1-Runtime ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/STL-Test1-Runtime )


# Platform-specific compiler options
if(MSVC_VERSION)
	set_target_properties (STL-Test1-Runtime PROPERTIES COMPILE_FLAGS ${MSVC_COMPILER_FLAGS})
	set_target_properties (STL-Test1-Runtime PROPERTIES LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)