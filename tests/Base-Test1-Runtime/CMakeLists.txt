###############################################################################
#
# CMakeLists.txt file for the GridRipper Base-Test1-Runtime test
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 23/07/2014
#
###############################################################################

# The project name and version, release numbers
project (Base-Test1-Runtime)


# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
include_directories (${PROJECT_SOURCE_DIR}/inc/)
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (Base_Test1_Runtime_BUILD)
set (Base_Test1_Runtime_HDRS)
set (Base_Test1_Runtime_SRCS)

set (Base_Test1_Runtime_HDRS ${Base_Test1_Runtime_HDRS} ${PROJECT_SOURCE_DIR}/inc/Base-Test1-Runtime.hpp)

set (Base_Test1_Runtime_SRCS ${Base_Test1_Runtime_SRCS} ${PROJECT_SOURCE_DIR}/src/Base-Test1-Runtime.cpp)

set (Base_Test1_Runtime_BUILD ${Base_Test1_Runtime_BUILD} ${Base_Test1_Runtime_HDRS} ${Base_Test1_Runtime_SRCS})


# Build and link
add_executable( Base-Test1-Runtime ${Base_Test1_Runtime_BUILD})
target_link_libraries( Base-Test1-Runtime Gripper ${Gripper_ADDITIONAL_LIBS})


# Create source file groups for IDEs
source_group ("Headers" FILES ${Base_Test1_Runtime_HDRS})
source_group ("Sources" FILES ${Base_Test1_Runtime_SRCS})

# Create project groups for IDEs
set_target_properties (Base-Test1-Runtime PROPERTIES FOLDER "Tests//Base")


# Add CTest entry
add_test( Base-Test1-Runtime ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Base-Test1-Runtime )


# Platform-specific compiler options
if(MSVC_VERSION)
#set_target_properties (Base-Test1-Runtime PROPERTIES
#                       COMPILE_FLAGS ${MSVC_COMPILER_FLAGS}
#                       LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)