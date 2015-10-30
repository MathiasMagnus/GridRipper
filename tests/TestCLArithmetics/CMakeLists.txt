###############################################################################
#
# CMakeLists.txt file for the GridRipper TestCLArithmetics test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2014/07/02
#
###############################################################################

# The project name and version, release numbers
project (TestCLArithmetics)


# Append our project's include directory to the "#include <>" paths
include_directories ("${PROJECT_SOURCE_DIR}/inc/")
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (TestCLArithmetics_BUILD)
set (TestCLArithmetics_HDRS)
set (TestCLArithmetics_SRCS)

set (TestCLArithmetics_HDRS ${TestCLArithmetics_HDRS} ${PROJECT_SOURCE_DIR}/inc/TestCLArithmetics.hpp)

set (TestCLArithmetics_SRCS ${TestCLArithmetics_SRCS} ${PROJECT_SOURCE_DIR}/src/TestCLArithmetics.cpp)

set (TestCLArithmetics_BUILD ${TestCLArithmetics_BUILD} ${TestCLArithmetics_HDRS} ${TestCLArithmetics_SRCS})

source_group ("Headers" FILES ${TestCLArithmetics_HDRS})
source_group ("Sources" FILES ${TestCLArithmetics_SRCS})

add_executable( TestCLArithmetics ${TestCLArithmetics_BUILD})
target_link_libraries( TestCLArithmetics Gripper ${Gripper_ADDITIONAL_LIBS})

add_test( TestCLArithmetics "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TestCLArithmetics -v")

if(MSVC_VERSION)
	set_target_properties (TestCLArithmetics PROPERTIES COMPILE_FLAGS ${MSVC_COMPILER_FLAGS})
	set_target_properties (TestCLArithmetics PROPERTIES LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)