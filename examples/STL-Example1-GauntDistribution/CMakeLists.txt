###############################################################################
#
# CMakeLists.txt file for the GridRipper STL-Example1-GauntDistribution example
# Hacker: Máté Ferenc Nagy-Egri
# Last modified: 2014/08/01
#
###############################################################################

# The project name and version, release numbers
project (STL-Example1-GauntDistribution)


# Append our project's include directory to the "#include <>" paths and inherit Gripper's own
include_directories (${PROJECT_SOURCE_DIR}/inc/)
include_directories (${Gripper_ADDITIONAL_INCLUDE_DIRS})


# Adding source code files according to configuration
set (STL_Example1_GauntDistribution_BUILD)
set (STL_Example1_GauntDistribution_HDRS)
set (STL_Example1_GauntDistribution_SRCS)

set (STL_Example1_GauntDistribution_HDRS ${STL_Example1_GauntDistribution_HDRS} ${PROJECT_SOURCE_DIR}/inc/STL_Example1_GauntDistribution.hpp)

set (STL_Example1_GauntDistribution_SRCS ${STL_Example1_GauntDistribution_SRCS} ${PROJECT_SOURCE_DIR}/src/STL_Example1_GauntDistribution.cpp)

set (STL_Example1_GauntDistribution_BUILD ${STL_Example1_GauntDistribution_BUILD} ${STL_Example1_GauntDistribution_HDRS} ${STL_Example1_GauntDistribution_SRCS})


# Adding library target for build
add_executable (STL-Example1-GauntDistribution ${STL_Example1_GauntDistribution_BUILD})
target_link_libraries (STL-Example1-GauntDistribution Gripper ${Gripper_ADDITIONAL_LIBS})


# Create project groups for IDEs
set_target_properties (STL-Example1-GauntDistribution PROPERTIES FOLDER "Examples//STL")

# Create source file groups for IDEs
source_group ("Headers" FILES ${STL_Example1_GauntDistribution_HDRS})
source_group ("Sources" FILES ${STL_Example1_GauntDistribution_SRCS})


# Add platform specific compile options
if(MSVC_VERSION)
#set_target_properties (STL_Example1_GauntDistribution PROPERTIES
#                       COMPILE_FLAGS ${MSVC_COMPILER_FLAGS}
#                       LINK_FLAGS ${MSVC_LINKER_FLAGS})
endif(MSVC_VERSION)