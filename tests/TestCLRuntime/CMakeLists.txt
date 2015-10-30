###############################################################################
#
# CMakeLists.txt file for the GridRipper TestCLRuntime test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2014/07/02
#
###############################################################################

# The project name and version, release numbers
project (TestCLRuntime)


# Append our project's include directory to the "#include <>" paths
include_directories ("${PROJECT_SOURCE_DIR}/inc/")
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (TestCLRuntime_BUILD)
set (TestCLRuntime_HDRS)
set (TestCLRuntime_SRCS)

set (TestCLRuntime_HDRS ${TestCLRuntime_HDRS} ${PROJECT_SOURCE_DIR}/inc/TestCLRuntime.hpp)

set (TestCLRuntime_SRCS ${TestCLRuntime_SRCS} ${PROJECT_SOURCE_DIR}/src/TestCLRuntime.cpp)

set (TestCLRuntime_BUILD ${TestCLRuntime_BUILD} ${TestCLRuntime_HDRS} ${TestCLRuntime_SRCS})

source_group ("Headers" FILES ${TestCLRuntime_HDRS})
source_group ("Sources" FILES ${TestCLRuntime_SRCS})

add_executable( TestCLRuntime ${TestCLRuntime_BUILD})
target_link_libraries( TestCLRuntime Gripper ${Gripper_ADDITIONAL_LIBS})

add_test( TestCLRuntime ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TestCLRuntime )

if(MSVC_VERSION)
	set_target_properties (TestCLRuntime PROPERTIES COMPILE_FLAGS ${MSVC_COMPILER_FLAGS})
	set_target_properties (TestCLRuntime PROPERTIES LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)